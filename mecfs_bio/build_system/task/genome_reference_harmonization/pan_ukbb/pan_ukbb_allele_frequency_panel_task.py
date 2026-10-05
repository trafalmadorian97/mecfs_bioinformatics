"""
Allele-frequency panel from the Pan-UK Biobank variant manifest, for genome-reference
harmonization.

The manifest (full_variant_qc_metrics.txt.bgz) is a bgzipped TSV of UK Biobank imputed
variants on GRCh37, already filtered to imputation INFO > 0.8, with the alternate allele
frequency af_{POP} for the groups AFR, AMR, CSA, EAS, EUR and MID. This Task streams it once
and writes CHR, POS, REF, ALT and one AF_ukb_{pop} column per configured group, through the
shared panel batch checks.

Every manifest row is kept except those whose ref disagrees with the FASTA. The manifest
holds a known set of such rows: ref and alt swapped, inherited from UK Biobank's imputed
BGEN allele order. expected_ref_mismatches pins their number, so any change in the file or in
the parsing fails the build. Allele numbers are not stored: in the manifest they are
constant within a contig and carry no per-site information.

The manifest is read through the asset's read_spec with polars, which detects the
compression, and streamed in batches of batch_rows rows. collect_batches is marked unstable
in polars; the Task tests pin that every row of a multi-block bgzipped file arrives.

See experiments/claude/design_specs/2026-10-02-gnomad-allele-frequency-panel-design.md.
"""

from collections.abc import Iterator, Sequence
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    gwaslab_code_to_contig_name,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    load_fasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    DEFAULT_BATCH_ROWS,
    SOURCE_CONTIG_COL,
    PanelBatchContext,
    PanelWriteSummary,
    write_checked_panel,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanUkbbGroup,
    panel_af_col,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL

logger = structlog.get_logger()

PAN_UKBB_BUILD: GenomeBuild = "19"
PAN_UKBB_PANEL_FILENAME = "pan_ukbb_allele_frequencies"
_MANIFEST_CONTIG_COL = "chrom"
_MANIFEST_POS_COL = "pos"
_MANIFEST_REF_COL = "ref"
_MANIFEST_ALT_COL = "alt"
# How to read the manifest. chrom must be read as a string: schema inference over the
# leading rows would otherwise type it as an integer and fail at "X".
PAN_UKBB_MANIFEST_READ_SPEC = DataFrameReadSpec(
    DataFrameTextFormat(
        separator="\t",
        null_values=["NA"],
        schema_overrides={_MANIFEST_CONTIG_COL: pl.String()},
    )
)


def manifest_af_col(group: PanUkbbGroup) -> str:
    """The manifest's frequency column for a group: ukb_eur -> af_EUR."""
    return "af_" + group.removeprefix("ukb_").upper()


@frozen(slots=True)
class PanUkbbAlleleFrequencyPanelTask(Task):
    meta: HarmonizableReferenceTableMeta
    manifest_task: Task
    fasta_task: Task
    groups: tuple[PanUkbbGroup, ...]
    chromosomes: tuple[int, ...]
    expected_ref_mismatches: int
    batch_rows: int = DEFAULT_BATCH_ROWS

    def __attrs_post_init__(self):
        assert self.groups, "at least one Pan-UKBB group is required"
        assert len(set(self.groups)) == len(self.groups), (
            f"duplicate groups {self.groups}"
        )
        assert self.chromosomes, "at least one chromosome is required"
        assert self.expected_ref_mismatches >= 0
        assert self.batch_rows > 0

    @property
    def deps(self) -> list[Task]:
        return [self.manifest_task, self.fasta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        manifest = scan_manifest(
            fetch(self.manifest_task.asset_id), self.manifest_task.meta
        )
        fasta = load_fasta(fetch, self.fasta_task)
        assert_manifest_columns(manifest, self.groups)
        out_path = scratch_dir / (PAN_UKBB_PANEL_FILENAME + ".parquet")
        summary = write_checked_panel(
            manifest_batches(manifest, self.groups, self.batch_rows),
            out_path,
            PanelBatchContext(
                fasta=fasta, contig_codes=manifest_contig_codes(self.chromosomes)
            ),
            byte_stream_split_columns=[],
        )
        assert_pan_ukbb_summary(summary, self.expected_ref_mismatches, self.chromosomes)
        logger.info(
            "Pan-UKBB allele-frequency panel written",
            rows_in=summary.counts.rows_in,
            rows_written=summary.rows_written,
            ref_mismatches_dropped=summary.counts.ref_mismatch,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        manifest_task: Task,
        fasta_task: Task,
        groups: Sequence[PanUkbbGroup],
        chromosomes: Sequence[int],
        expected_ref_mismatches: int,
    ) -> "PanUkbbAlleleFrequencyPanelTask":
        fasta_meta = fasta_task.meta
        assert isinstance(fasta_meta, FASTAMeta), (
            f"fasta_task must carry FASTAMeta, got {type(fasta_meta).__name__}"
        )
        assert fasta_meta.build == PAN_UKBB_BUILD, (
            f"the Pan-UKBB manifest is GRCh37; got a build {fasta_meta.build} FASTA"
        )
        source_meta = manifest_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta manifest, got {type(source_meta).__name__}"
        )
        assert source_meta.read_spec is not None, (
            "the manifest meta needs a read_spec (PAN_UKBB_MANIFEST_READ_SPEC)"
        )
        return cls(
            meta=HarmonizableReferenceTableMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename=PAN_UKBB_PANEL_FILENAME,
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=PAN_UKBB_BUILD,
                    ref_allele_col=PANEL_REF_COL,
                    pos_col=GWASLAB_POS_COL,
                ),
                allele_frequency_columns=PanelAlleleFrequencyColumns.prefixed(groups),
            ),
            manifest_task=manifest_task,
            fasta_task=fasta_task,
            groups=tuple(groups),
            chromosomes=tuple(chromosomes),
            expected_ref_mismatches=expected_ref_mismatches,
        )


def manifest_contig_codes(chromosomes: Sequence[int]) -> dict[str, int]:
    """Manifest contig names ("1".."22", "X") for the expected gwaslab codes."""
    return {gwaslab_code_to_contig_name(code): code for code in chromosomes}


def scan_manifest(asset: Asset, meta: Meta) -> pl.LazyFrame:
    """The manifest as a polars LazyFrame, read through the asset's read_spec."""
    native = scan_dataframe_asset(asset, meta).to_native()
    assert isinstance(native, pl.LazyFrame), (
        f"expected a polars LazyFrame for the manifest, got {type(native).__name__}"
    )
    return native


def assert_manifest_columns(
    manifest: pl.LazyFrame, groups: Sequence[PanUkbbGroup]
) -> None:
    present = manifest.collect_schema().names()
    required = [
        _MANIFEST_CONTIG_COL,
        _MANIFEST_POS_COL,
        _MANIFEST_REF_COL,
        _MANIFEST_ALT_COL,
        *[manifest_af_col(group) for group in groups],
    ]
    missing = [column for column in required if column not in present]
    assert not missing, f"the Pan-UKBB manifest lacks columns {missing}"


def manifest_batches(
    manifest: pl.LazyFrame, groups: Sequence[PanUkbbGroup], batch_rows: int
) -> Iterator[pl.DataFrame]:
    """Stream the manifest's needed columns, renamed to the panel's, in batches."""
    af_columns = {manifest_af_col(group): panel_af_col(group) for group in groups}
    query = manifest.select(
        pl.col(_MANIFEST_CONTIG_COL).cast(pl.String).alias(SOURCE_CONTIG_COL),
        pl.col(_MANIFEST_POS_COL).alias(GWASLAB_POS_COL),
        pl.col(_MANIFEST_REF_COL).alias(PANEL_REF_COL),
        pl.col(_MANIFEST_ALT_COL).alias(PANEL_ALT_COL),
        *[
            pl.col(source).cast(pl.Float32).alias(target)
            for source, target in af_columns.items()
        ],
    )
    for batch in query.collect_batches(chunk_size=batch_rows):
        null_counts = batch.select(
            [pl.col(column).null_count() for column in af_columns.values()]
        ).row(0, named=True)
        assert not any(null_counts.values()), (
            f"null allele frequencies in the Pan-UKBB manifest: {null_counts}"
        )
        yield batch


def assert_pan_ukbb_summary(
    summary: PanelWriteSummary,
    expected_ref_mismatches: int,
    chromosomes: Sequence[int],
) -> None:
    assert summary.counts.ref_mismatch == expected_ref_mismatches, (
        f"expected {expected_ref_mismatches} manifest rows whose ref differs from the "
        f"FASTA, found {summary.counts.ref_mismatch}"
    )
    assert summary.counts.fasta_ambiguous == 0, (
        f"{summary.counts.fasta_ambiguous} manifest rows lie over N or IUPAC FASTA bases"
    )
    missing = sorted(set(chromosomes) - summary.chromosomes)
    assert not missing, f"no Pan-UKBB rows on chromosomes {missing}"
    assert summary.rows_written > 0, "the Pan-UKBB panel is empty"
