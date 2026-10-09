"""
Invariants shared by the allele-frequency panel extraction Tasks, applied batch by batch.

A panel source (a gnomAD sites VCF, the Pan-UKBB variant manifest) is streamed as polars
batches holding a source contig column, POS, REF, ALT and frequency columns. Each batch is
checked and converted here, and the surviving rows are written to one parquet file.

Fatal checks:
- contig names must be among the expected ones; they become an Int32 gwaslab CHR code;
- REF and ALT must be non-null and contain only A, C, G and T;
- (CHR, POS) must never decrease, and no (CHR, POS, REF, ALT) key may repeat. The keys at
  the last position of each batch are carried into the next, so a pair split across a batch
  boundary is still caught;
- every REF span must lie inside its FASTA contig (a span past the end means a wrong build or
  contig naming).

REF is then compared with the FASTA. Records whose reference span holds N or another IUPAC
code are dropped and counted as fasta_ambiguous: their orientation is unverifiable, and the
harmonizer's own FASTA classification would reject them anyway. Records whose REF differs
from a pure A/C/G/T span are dropped and counted as ref_mismatch. The calling Task decides
how many of each it accepts.
"""

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    open_parquet_writer,
    parquet_encoding,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    DEFAULT_MAX_GATHER_BYTES,
    IndexedFasta,
    reference_is_acgt,
    reference_matches,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

SOURCE_CONTIG_COL = "contig"
PANEL_KEY_COLUMNS = [GWASLAB_CHROM_COL, GWASLAB_POS_COL, PANEL_REF_COL, PANEL_ALT_COL]
# Rows per streamed batch; each batch becomes one parquet row group. On the full Pan-UKBB
# manifest, 1,000,000-row batches peaked at about 3 GB RSS.
DEFAULT_BATCH_ROWS = 1_000_000
_ACGT_ALLELE = r"^[ACGT]+$"
_FASTA_ACGT_COL = "_fasta_acgt"
_REF_MATCHES_COL = "_ref_matches"
_CHROM_STEP_COL = "_chrom_step"
_POS_STEP_COL = "_pos_step"


@frozen(slots=True)
class PanelBatchCounts:
    rows_in: int
    fasta_ambiguous: int
    ref_mismatch: int

    @classmethod
    def zero(cls) -> "PanelBatchCounts":
        return cls(rows_in=0, fasta_ambiguous=0, ref_mismatch=0)

    def __add__(self, other: "PanelBatchCounts") -> "PanelBatchCounts":
        return PanelBatchCounts(
            rows_in=self.rows_in + other.rows_in,
            fasta_ambiguous=self.fasta_ambiguous + other.fasta_ambiguous,
            ref_mismatch=self.ref_mismatch + other.ref_mismatch,
        )


@frozen(slots=True)
class PanelBatchContext:
    """What every batch is checked against: the FASTA and the allowed source contigs."""

    fasta: IndexedFasta
    contig_codes: Mapping[str, int]
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES


@frozen(slots=True)
class CheckedPanelBatch:
    table: pl.DataFrame
    counts: PanelBatchCounts
    carry: pl.DataFrame


@frozen(slots=True)
class PanelWriteSummary:
    rows_written: int
    counts: PanelBatchCounts
    chromosomes: frozenset[int]


def fetch_file_path(fetch: Fetch, task: Task) -> Path:
    """Fetch a Task's FileAsset and return its path."""
    asset = fetch(task.asset_id)
    assert isinstance(asset, FileAsset), f"expected {task.asset_id} to be a FileAsset"
    return asset.path


def write_checked_panel(
    batches: Iterable[pl.DataFrame],
    out_path: Path,
    context: PanelBatchContext,
    byte_stream_split_columns: Sequence[str],
) -> PanelWriteSummary:
    """Check every batch, write the kept rows to out_path as one parquet, and summarize."""
    counts = PanelBatchCounts.zero()
    carry = _empty_keys()
    rows_written = 0
    chromosomes: set[int] = set()
    writer: pyarrow.parquet.ParquetWriter | None = None
    try:
        for batch in batches:
            checked = check_panel_batch(batch, context, carry)
            carry = checked.carry
            counts = counts + checked.counts
            table = checked.table.to_arrow()
            if writer is None:
                writer = open_parquet_writer(
                    out_path,
                    table.schema,
                    parquet_encoding(
                        table.schema.names,
                        compression="zstd",
                        compression_level=None,
                        byte_stream_split_columns=byte_stream_split_columns,
                    ),
                )
            # An empty row group gets zero page offsets, which polars refuses to read.
            if table.num_rows > 0:
                writer.write_table(table)
            rows_written += checked.table.height
            chromosomes.update(checked.table[GWASLAB_CHROM_COL].unique().to_list())
    finally:
        if writer is not None:
            writer.close()
    assert writer is not None, "the panel source produced no batches"
    return PanelWriteSummary(
        rows_written=rows_written, counts=counts, chromosomes=frozenset(chromosomes)
    )


def check_panel_batch(
    batch: pl.DataFrame, context: PanelBatchContext, carry: pl.DataFrame
) -> CheckedPanelBatch:
    """Apply the module's checks to one batch; carry holds the previous batch's last keys."""
    converted = _with_chromosome_codes(batch, context.contig_codes)
    _assert_acgt_alleles(converted)
    keys = pl.concat([carry, converted.select(PANEL_KEY_COLUMNS)])
    _assert_sorted_and_unique(keys)
    classified = _classify_against_fasta(converted, context)
    ambiguous = classified.select((~pl.col(_FASTA_ACGT_COL)).sum()).item()
    mismatched = classified.select(
        (pl.col(_FASTA_ACGT_COL) & ~pl.col(_REF_MATCHES_COL)).sum()
    ).item()
    return CheckedPanelBatch(
        table=classified.filter(
            pl.col(_FASTA_ACGT_COL) & pl.col(_REF_MATCHES_COL)
        ).drop(_FASTA_ACGT_COL, _REF_MATCHES_COL),
        counts=PanelBatchCounts(
            rows_in=converted.height,
            fasta_ambiguous=int(ambiguous),
            ref_mismatch=int(mismatched),
        ),
        carry=_keys_at_last_position(keys, carry),
    )


def _empty_keys() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
        }
    )


def _with_chromosome_codes(
    batch: pl.DataFrame, contig_codes: Mapping[str, int]
) -> pl.DataFrame:
    contigs = set(batch[SOURCE_CONTIG_COL].unique().to_list())
    unexpected = contigs - set(contig_codes)
    assert not unexpected, (
        f"unexpected contigs {sorted(map(str, unexpected))}; expected {sorted(contig_codes)}"
    )
    value_columns = [
        column
        for column in batch.columns
        if column
        not in (SOURCE_CONTIG_COL, GWASLAB_POS_COL, PANEL_REF_COL, PANEL_ALT_COL)
    ]
    return batch.select(
        pl.col(SOURCE_CONTIG_COL)
        .replace_strict(dict(contig_codes), return_dtype=pl.Int32)
        .alias(GWASLAB_CHROM_COL),
        pl.col(GWASLAB_POS_COL).cast(pl.Int32),
        PANEL_REF_COL,
        PANEL_ALT_COL,
        *value_columns,
    )


def _assert_acgt_alleles(batch: pl.DataFrame) -> None:
    for column in (PANEL_REF_COL, PANEL_ALT_COL):
        bad = batch.filter(
            pl.col(column).is_null() | ~pl.col(column).str.contains(_ACGT_ALLELE)
        )
        assert bad.height == 0, (
            f"{bad.height} records with a null or non-ACGT {column}, e.g. "
            f"{bad.select(PANEL_KEY_COLUMNS).head(3).rows()}"
        )


def _assert_sorted_and_unique(keys: pl.DataFrame) -> None:
    steps = keys.with_columns(
        pl.col(GWASLAB_CHROM_COL).diff().alias(_CHROM_STEP_COL),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64).diff().alias(_POS_STEP_COL),
    )
    out_of_order = steps.filter(
        (pl.col(_CHROM_STEP_COL) < 0)
        | ((pl.col(_CHROM_STEP_COL) == 0) & (pl.col(_POS_STEP_COL) < 0))
    )
    assert out_of_order.height == 0, (
        "records out of (CHR, POS) order at "
        f"{out_of_order.select(GWASLAB_CHROM_COL, GWASLAB_POS_COL).head(3).rows()}"
    )
    duplicated = keys.filter(pl.struct(PANEL_KEY_COLUMNS).is_duplicated())
    assert duplicated.height == 0, (
        f"duplicate panel keys, e.g. {duplicated.unique().head(3).rows()}"
    )


def _keys_at_last_position(keys: pl.DataFrame, carry: pl.DataFrame) -> pl.DataFrame:
    if keys.height == 0:
        return carry
    last = keys.row(-1, named=True)
    return keys.filter(
        (pl.col(GWASLAB_CHROM_COL) == last[GWASLAB_CHROM_COL])
        & (pl.col(GWASLAB_POS_COL) == last[GWASLAB_POS_COL])
    )


def _classify_against_fasta(
    batch: pl.DataFrame, context: PanelBatchContext
) -> pl.DataFrame:
    if batch.height == 0:
        return batch.with_columns(
            pl.lit(True).alias(_FASTA_ACGT_COL), pl.lit(True).alias(_REF_MATCHES_COL)
        )
    return pl.concat(
        [
            _classify_chromosome(rows, int(code), context)
            for (code,), rows in batch.partition_by(
                GWASLAB_CHROM_COL, as_dict=True, maintain_order=True
            ).items()
        ]
    )


def _classify_chromosome(
    rows: pl.DataFrame, code: int, context: PanelBatchContext
) -> pl.DataFrame:
    fasta = context.fasta
    assert code in fasta.entries, f"chromosome {code} is not in {fasta.fasta_path}"
    positions = rows[GWASLAB_POS_COL].to_numpy()
    lengths = rows[PANEL_REF_COL].str.len_bytes().to_numpy()
    span_end = positions.astype(np.int64) + lengths - 1
    beyond = (positions < 1) | (span_end > fasta.entries[code].length)
    assert not beyond.any(), (
        f"chromosome {code}: {int(beyond.sum())} REF spans lie outside the FASTA contig "
        "(wrong build or contig naming?)"
    )
    return rows.with_columns(
        pl.Series(
            _FASTA_ACGT_COL,
            reference_is_acgt(
                fasta, code, positions, lengths, context.max_gather_bytes
            ),
        ),
        pl.Series(
            _REF_MATCHES_COL,
            reference_matches(
                fasta, code, positions, rows[PANEL_REF_COL], context.max_gather_bytes
            ),
        ),
    )
