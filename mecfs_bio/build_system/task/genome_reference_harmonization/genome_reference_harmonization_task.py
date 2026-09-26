"""
Genome-reference harmonization of GWAS summary statistics.

Orients every variant so that NEA is the plus-strand reference allele of a genome
FASTA, flips allele-dependent statistics to match, and resolves the variants the genome
alone cannot orient.

- **Trust.** A table whose checkable SNVs and indels are 100% reference-consistent,
  and with a low number of "suspicious indels" is trusted: its palindromic variants
  and its indels with both alleles on the genome keep their source orientation.
- **Untrusted tables.** Palindromic SNV strands and ambiguous indels are resolved against
  a reference panel's allele frequencies, and are dropped whenever the evidence is not
  decisive.

This replaces gwaslab harmonization. See the design spec in
experiments/claude/design_specs/2026-09-14-genome-reference-harmonization-design.md.

Memory is bounded by one chromosome. Pass 1 reads the key columns and, when present, EAF
per chromosome to decide trust, consulting the panel at ambiguous-indel positions for the
suspicious-indel gate. Pass 2 resolves one chromosome at a time into parquet parts, which
are then concatenated with a streaming sink.
"""

from collections.abc import Sequence
from pathlib import Path

import attrs
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_BOTH,
    classify_alleles,
    prepare_alleles,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    DROPPED_COLUMNS,
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
    resolve_chromosome,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    SuspiciousIndelCounts,
    TrustCounts,
    TrustEvidence,
    count_suspicious_indels,
    count_trust_evidence,
    decide_trust,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

logger = structlog.get_logger()

HARMONIZED_FILENAME = "harmonized.parquet"
# A pandas index serialized into a parquet by pyarrow; never a statistic, so dropped rather
# than failing the unregistered-column check.
_PANDAS_INDEX_COLUMN = "__index_level_0__"
_KEPT_LABEL = "kept"  # log label for rows without a drop reason
_KEY_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


def scan_sumstats_as_polars(
    asset: Asset, meta: Meta, pipe: DataProcessingPipe
) -> pl.LazyFrame:
    """Scan the sumstats asset, apply the pipe, and insist the result is still polars-backed."""
    native = pipe.process(scan_dataframe_asset(asset, meta)).to_native()
    assert isinstance(native, pl.LazyFrame), (
        "genome-reference harmonization streams polars LazyFrames, but the pipe produced "
        f"{type(native).__name__}"
    )
    if _PANDAS_INDEX_COLUMN in native.collect_schema().names():
        native = native.drop(_PANDAS_INDEX_COLUMN)
    return native


def chromosomes_to_harmonize(
    sumstats: pl.LazyFrame,
    fasta: IndexedFasta,
    options: GenomeReferenceHarmonizationOptions,
) -> list[int]:
    null_counts = (
        sumstats.select([pl.col(column).null_count() for column in _KEY_COLUMNS])
        .collect(engine="streaming")
        .row(0, named=True)
    )
    assert all(count == 0 for count in null_counts.values()), (
        f"null values in key columns: {null_counts}"
    )
    present = sorted(
        int(chrom)
        for chrom in sumstats.select(pl.col(GWASLAB_CHROM_COL).unique())
        .collect(engine="streaming")[GWASLAB_CHROM_COL]
        .to_list()
    )
    excluded = [chrom for chrom in present if chrom in options.excluded_chromosomes]
    if excluded:
        logger.info("dropping rows on excluded chromosomes", chromosomes=excluded)
    kept = [chrom for chrom in present if chrom not in options.excluded_chromosomes]
    missing = [chrom for chrom in kept if chrom not in fasta.entries]
    assert not missing, f"chromosomes {missing} are not in the FASTA {fasta.fasta_path}"
    return kept


@frozen(slots=True)
class ParquetPanelLoader:
    """Reads one chromosome's panel rows at requested positions from the panel parquet."""

    panel_path: Path
    chrom: int

    def __call__(self, positions: pl.Series) -> pl.DataFrame:
        wanted = positions.cast(pl.Int32).unique().to_frame(GWASLAB_POS_COL).lazy()
        return (
            pl.scan_parquet(self.panel_path)
            .filter(pl.col(GWASLAB_CHROM_COL) == self.chrom)
            .join(wanted, on=GWASLAB_POS_COL, how="semi")
            .select(
                pl.col(GWASLAB_POS_COL).cast(pl.Int64),
                PANEL_REF_COL,
                PANEL_ALT_COL,
                PANEL_AF_COL,
            )
            .collect()
        )


def _ambiguous_positions(classified: pl.DataFrame) -> pl.Series:
    """Positions of ambiguous (BOTH) indels, whose orientation the suspicion test checks."""
    return classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH)[
        GWASLAB_POS_COL
    ]


def count_trust_evidence_genome_wide(
    sumstats: pl.LazyFrame,
    chromosomes: Sequence[int],
    fasta: IndexedFasta,
    panel_path: Path,
    options: GenomeReferenceHarmonizationOptions,
) -> TrustEvidence:
    names = sumstats.collect_schema().names()
    eaf_present = GWASLAB_EFFECT_ALLELE_FREQ_COL in names
    columns = [
        column
        for column in [*_KEY_COLUMNS, GWASLAB_EFFECT_ALLELE_FREQ_COL]
        if column in names
    ]
    counts = TrustCounts.zero()
    suspicious = SuspiciousIndelCounts.zero()
    for chrom in chromosomes:
        frame = (
            sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
            .select(columns)
            .collect(engine="streaming")
        )
        valid = prepare_alleles(frame).filter(valid_alleles_expr())
        classified = classify_alleles(
            valid, fasta=fasta, chrom=chrom, max_gather_bytes=options.max_gather_bytes
        )
        counts = counts + count_trust_evidence(classified)
        if eaf_present:
            panel = ParquetPanelLoader(panel_path=panel_path, chrom=chrom)(
                _ambiguous_positions(classified)
            )
            suspicious = suspicious + count_suspicious_indels(
                classified, panel, options
            )
    return TrustEvidence(counts=counts, suspicious=suspicious, eaf_present=eaf_present)


def _variant_key_is_duplicated() -> pl.Expr:
    """Whether each row shares its within-chromosome (POS, EA, NEA) key with another row."""
    return pl.struct(
        GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
    ).is_duplicated()


def resolve_chromosome_rows(
    sumstats: pl.LazyFrame, context: ChromosomeContext, panel_path: Path
) -> pl.DataFrame:
    """All rows of one chromosome, resolved, with DROP_REASON_COL (used by experiments too)."""
    rows = sumstats.filter(pl.col(GWASLAB_CHROM_COL) == context.chrom).collect(
        engine="streaming"
    )
    # Defensive: the input must already be unique per variant key. Mirrored pairs (distinct
    # keys) are fine and preserved; genuine duplicate keys are a source-data error and are
    # removed upstream (transform_gwaslab_sumstats), so reaching here means that failed.
    input_duplicated = rows.select(_variant_key_is_duplicated().any()).item()
    assert not input_duplicated, (
        f"chromosome {context.chrom}: input summary statistics have duplicate "
        "(CHR, POS, EA, NEA) keys; harmonization requires unique input variants"
    )
    return resolve_chromosome(
        rows,
        context,
        load_panel=ParquetPanelLoader(panel_path=panel_path, chrom=context.chrom),
    )


def _write_chromosome_part(
    sumstats: pl.LazyFrame,
    context: ChromosomeContext,
    panel_path: Path,
    parts_dir: Path,
) -> Path:
    resolved = resolve_chromosome_rows(sumstats, context, panel_path)
    rows_by_reason = {
        (_KEPT_LABEL if reason is None else reason): count
        for reason, count in resolved.group_by(DROP_REASON_COL).len().rows()
    }
    logger.info(
        "genome-reference harmonization chromosome summary",
        chromosome=context.chrom,
        rows_by_drop_reason=rows_by_reason,
    )
    kept = (
        resolved.filter(pl.col(DROP_REASON_COL).is_null())
        .drop(
            DROP_REASON_COL,
            *[column for column in DROPPED_COLUMNS if column in resolved.columns],
        )
        .sort(GWASLAB_POS_COL)
    )
    # Orientation can map two distinct input variants onto the same (CHR, POS, EA, NEA) key
    # (a liftover many-to-one collision, or a mirrored pair that resolves to one orientation).
    # Such a key is ambiguous for any allele-keyed join, and we cannot tell which variant is
    # correct, so every colliding row is dropped. Mirrored pairs that stay in distinct
    # orientations keep distinct keys and are not affected.
    collision = kept.select(_variant_key_is_duplicated()).to_series()
    collision_count = int(collision.sum())
    if collision_count:
        logger.warning(
            "dropping variants whose oriented key collides after harmonization",
            chromosome=context.chrom,
            count=collision_count,
        )
        kept = kept.filter(~collision)
    part_path = parts_dir / f"chr{context.chrom}.parquet"
    kept.write_parquet(part_path)
    return part_path


def load_fasta(fetch: Fetch, fasta_task: Task) -> IndexedFasta:
    """Fetch the FASTA directory asset and open it as an IndexedFasta."""
    fasta_asset = fetch(fasta_task.asset_id)
    assert isinstance(fasta_asset, DirectoryAsset), (
        f"expected {fasta_task.asset_id} to be a DirectoryAsset"
    )
    return IndexedFasta.open(fasta_asset.path)


def load_panel_path(fetch: Fetch, panel_task: Task) -> Path:
    """Fetch the reference panel file asset and return its path."""
    panel_asset = fetch(panel_task.asset_id)
    assert isinstance(panel_asset, FileAsset), (
        f"expected {panel_task.asset_id} to be a FileAsset"
    )
    return panel_asset.path


def resolve_harmonized_build(fasta_task: Task, panel_task: Task) -> GenomeBuild:
    """Assert the FASTA and the allele-frequency panel are the same genome build and
    return it."""
    fasta_meta = fasta_task.meta
    panel_meta = panel_task.meta
    assert isinstance(fasta_meta, FASTAMeta), (
        f"fasta_task must carry FASTAMeta, got {type(fasta_meta).__name__}"
    )
    assert (
        isinstance(panel_meta, HarmonizableReferenceTableMeta)
        and panel_meta.harmonization_info is not None
    ), "panel_task must carry HarmonizableReferenceTableMeta with harmonization_info"
    assert fasta_meta.build == panel_meta.harmonization_info.build, (
        f"fasta build {fasta_meta.build} != panel build {panel_meta.harmonization_info.build}"
    )
    return fasta_meta.build


@frozen(slots=True)
class GenomeReferenceHarmonizationTask(Task):
    meta: FilteredGWASDataMeta
    sumstats_task: Task
    fasta_task: Task
    panel_task: Task
    options: GenomeReferenceHarmonizationOptions
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def deps(self) -> list[Task]:
        return [self.sumstats_task, self.fasta_task, self.panel_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats = scan_sumstats_as_polars(
            fetch(self.sumstats_task.asset_id), self.sumstats_task.meta, self.pipe
        )
        fasta = load_fasta(fetch, self.fasta_task)
        panel_path = load_panel_path(fetch, self.panel_task)
        rules = resolve_column_rules(
            columns=sumstats.collect_schema().names(),
            extra=self.options.extra_column_rules,
        )
        chromosomes = chromosomes_to_harmonize(sumstats, fasta, self.options)
        evidence = count_trust_evidence_genome_wide(
            sumstats, chromosomes, fasta, panel_path, self.options
        )
        trusted = decide_trust(evidence, self.options)
        logger.info(
            "genome-reference harmonization trust decision",
            trusted=trusted,
            eaf_present=evidence.eaf_present,
            counts=attrs.asdict(evidence.counts),
            suspicious=attrs.asdict(evidence.suspicious),
        )
        parts_dir = scratch_dir / "parts"
        parts_dir.mkdir()
        part_paths = [
            _write_chromosome_part(
                sumstats,
                ChromosomeContext(
                    chrom=chrom,
                    fasta=fasta,
                    trusted=trusted,
                    rules=rules,
                    options=self.options,
                ),
                panel_path=panel_path,
                parts_dir=parts_dir,
            )
            for chrom in chromosomes
        ]
        assert part_paths, "no chromosomes to harmonize"
        out_path = scratch_dir / HARMONIZED_FILENAME
        pl.concat([pl.scan_parquet(path) for path in part_paths]).sink_parquet(out_path)
        n_rows = pl.scan_parquet(out_path).select(pl.len()).collect().item()
        assert n_rows > 0, "no variants survived genome-reference harmonization"
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sumstats_task: Task,
        fasta_task: Task,
        panel_task: Task,
        options: GenomeReferenceHarmonizationOptions = GenomeReferenceHarmonizationOptions(),
        pipe: DataProcessingPipe = IdentityPipe(),
    ) -> "GenomeReferenceHarmonizationTask":
        source_meta = sumstats_task.meta
        assert isinstance(source_meta, FilteredGWASDataMeta), (
            f"expected a FilteredGWASDataMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        build = resolve_harmonized_build(fasta_task, panel_task)
        return cls(
            meta=FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir="processed",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=build,
                    ref_allele_col=GWASLAB_NON_EFFECT_ALLELE_COL,
                    pos_col=GWASLAB_POS_COL,
                ),
            ),
            sumstats_task=sumstats_task,
            fasta_task=fasta_task,
            panel_task=panel_task,
            options=options,
            pipe=pipe,
        )
