"""
Task to combine GWAS with fixed-effects meta analysis
"""

from pathlib import Path, PurePath
from typing import Self, Sequence

import narwhals
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
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
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

logger = structlog.get_logger()

# Number of studies contributing to each output variant. Not a gwaslab column, so it
# lives here rather than in the gwaslab constants.
META_ANALYSIS_NUM_STUDIES_COL = "N_STUDIES"

# Columns that exist only between preparation and aggregation, never in the output.
_WEIGHT_COL = "_weight"
_WEIGHTED_BETA_COL = "_weighted_beta"
_SOURCE_SAMPLE_SIZE_COL = "_source_sample_size"

# Names of the per-source validation counts.
_NUM_ROWS = "_num_rows"
_NUM_DISTINCT_VARIANTS = "_num_distinct_variants"
_NUM_USABLE_SE = "_num_usable_se"
_NUM_USABLE_BETA = "_num_usable_beta"
_NUM_SELF_REFERENTIAL = "_num_self_referential"

# Every source is harmonized so that its non-effect allele is the reference allele, so a
# variant is identified exactly by these columns, with no allele reordering.
_VARIANT_KEY_COLS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_COL,
]


@frozen(slots=True)
class CaseControlSampleInfo:
    """
    Specifies the number of cases and controls in a case control study
    """

    cases: int
    controls: int

    def effective_sample_size(self) -> int:
        return int(4 / (1 / self.cases + 1 / self.controls))


SampleInfo = CaseControlSampleInfo  # add other types of sample info later


@frozen(slots=True)
class GwasSource:
    """
    Describes a source from which to draw GWAS data in order to perform a meta analysis

    The pipe is applied after reading and must not change allele orientation, since the
    meta analysis relies on the task's harmonization to the reference genome.
    """

    task: Task
    sample_info: SampleInfo
    pipe: DataProcessingPipe = IdentityPipe()


@frozen(slots=True)
class FixedEffectsMetaAnalysisTask(Task):
    """
    Task to perform a fixed effects meta analysis on non-overlapping GWAS of the same trait

    Every source must be harmonized to the same reference genome build, so that its
    non-effect allele is the reference allele (as produced by
    GenomeReferenceHarmonizationTask). This is checked from the sources' metadata at
    construction. A variant is then identified exactly by chromosome, position, reference
    allele and effect allele, and no allele flipping is performed. Identifying variants
    by their unordered allele pair instead would be unsafe: an insertion and the mirrored
    deletion at the same position (for example T/TCA and TCA/T) share an allele pair but
    are distinct variants.

    The variants present in the output dataframe are the union of the variants in the
    studies, so a variant genotyped by only some of the studies is kept and is combined
    from those studies alone.

    Each output variant carries a sample size column, GWASLAB_SAMPLE_SIZE_COLUMN, holding
    the sum of the effective sample sizes of just the studies that contributed to it, and
    a count of those studies in META_ANALYSIS_NUM_STUDIES_COL. Because that sample size
    varies per variant, downstream analyses should consume it with a PerVariantSampleSize
    rather than a scalar.

    Output rows are ordered by chromosome and position.

    For more information on fixed effects meta analysis, see:
    Chapter 22 of
    Balding, David J., Ida Moltke, and John Marioni, eds. Handbook of statistical genomics. John Wiley & Sons, 2019.
    """

    meta: Meta
    sources: Sequence[GwasSource]

    def __attrs_post_init__(self) -> None:
        assert len(self.sources) > 1
        common_harmonized_build(self.sources)

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.sources]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_dfs = [read_source(source, fetch) for source in self.sources]
        for source_index, df in enumerate(source_dfs):
            check_source(df, source_index=source_index)

        keep_rsid = any(GWASLAB_RSID_COL in df.collect_schema() for df in source_dfs)
        prepared_dfs = [
            prepare_source(df, sample_info=source.sample_info, keep_rsid=keep_rsid)
            for df, source in zip(source_dfs, self.sources)
        ]

        combined = combine_prepared_sources(prepared_dfs, keep_rsid=keep_rsid)
        out_path = scratch_dir / (self.asset_id + ".parquet")
        combined.sink_parquet(out_path)

        final = narwhals.scan_parquet(out_path, backend="polars")
        report_output_size(final)
        report_study_count_distribution(final)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, meta_analysis_name: str, sources: Sequence[GwasSource]
    ) -> Self:
        assert len(sources) > 1
        source_meta = sources[0].task.meta
        assert isinstance(source_meta, FilteredGWASDataMeta)
        for item in sources[1:]:
            assert isinstance(item.task.meta, FilteredGWASDataMeta)
            assert item.task.meta.trait == source_meta.trait
        build = common_harmonized_build(sources)
        meta = FilteredGWASDataMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=meta_analysis_name,
            sub_dir=PurePath("processed"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            genome_build=build,
            harmonization_info=HarmonizationInfo(
                build=build,
                ref_allele_col=GWASLAB_NON_EFFECT_ALLELE_COL,
                pos_col=GWASLAB_POS_COL,
            ),
        )
        return cls(
            meta=meta,
            sources=sources,
        )


def common_harmonized_build(sources: Sequence[GwasSource]) -> GenomeBuild:
    """The genome build that every source is harmonized to.

    Asserts that each source's metadata records harmonization with the non-effect allele
    as the reference allele, and that all sources share one build.
    """
    builds: set[GenomeBuild] = set()
    for source_index, source in enumerate(sources):
        source_meta = source.task.meta
        assert (
            isinstance(source_meta, FilteredGWASDataMeta)
            and source_meta.harmonization_info is not None
        ), (
            f"Source {source_index} ({source.task.asset_id}) must be a "
            "FilteredGWASDataMeta with harmonization_info; harmonize it with "
            "GenomeReferenceHarmonizationTask first."
        )
        info = source_meta.harmonization_info
        assert (
            info.ref_allele_col == GWASLAB_NON_EFFECT_ALLELE_COL
            and info.pos_col == GWASLAB_POS_COL
        ), (
            f"Source {source_index} ({source.task.asset_id}) is harmonized with "
            f"reference allele column {info.ref_allele_col} and position column "
            f"{info.pos_col}, but the meta analysis requires {GWASLAB_NON_EFFECT_ALLELE_COL} "
            f"and {GWASLAB_POS_COL}."
        )
        builds.add(info.build)
    assert len(builds) == 1, (
        f"Sources are harmonized to different genome builds: {sorted(builds)}"
    )
    return builds.pop()


def read_source(source: GwasSource, fetch: Fetch) -> narwhals.LazyFrame:
    asset = fetch(source.task.asset_id)
    return source.pipe.process(scan_dataframe_asset(asset, source.task.meta))


def check_source(df: narwhals.LazyFrame, source_index: int) -> None:
    """Fail fast on input that a fixed effects meta analysis cannot combine.

    Every check is gathered in one streaming pass, so validating a study costs one scan
    of it however many things are being checked.
    """
    counts = df.select(
        narwhals.len().alias(_NUM_ROWS),
        narwhals.concat_str(
            [narwhals.col(c).cast(narwhals.String()) for c in _VARIANT_KEY_COLS],
            separator=":",
        )
        .n_unique()
        .alias(_NUM_DISTINCT_VARIANTS),
        # A null or NaN fails the comparison too, so this one check covers zero,
        # negative and missing standard errors alike.
        (
            (narwhals.col(GWASLAB_SE_COL).cast(narwhals.Float64()) > 0)
            & narwhals.col(GWASLAB_SE_COL).cast(narwhals.Float64()).is_finite()
        )
        .sum()
        .alias(_NUM_USABLE_SE),
        narwhals.col(GWASLAB_BETA_COL)
        .cast(narwhals.Float64())
        .is_finite()
        .sum()
        .alias(_NUM_USABLE_BETA),
        (
            narwhals.col(GWASLAB_EFFECT_ALLELE_COL)
            == narwhals.col(GWASLAB_NON_EFFECT_ALLELE_COL)
        )
        .sum()
        .alias(_NUM_SELF_REFERENTIAL),
    ).collect(engine="streaming")

    num_rows = counts.item(0, _NUM_ROWS)
    assert counts.item(0, _NUM_DISTINCT_VARIANTS) == num_rows, (
        f"Source {source_index} states some variant more than once. Its effect size "
        "would be counted more than once by the meta analysis."
    )
    assert counts.item(0, _NUM_USABLE_SE) == num_rows, (
        f"Source {source_index} has standard errors that are not finite and positive. "
        "The meta analysis weights each study by the inverse of its squared standard "
        "error."
    )
    assert counts.item(0, _NUM_USABLE_BETA) == num_rows, (
        f"Source {source_index} has effect sizes that are not finite."
    )
    assert counts.item(0, _NUM_SELF_REFERENTIAL) == 0, (
        f"Source {source_index} has variants whose effect and non-effect alleles are "
        "the same."
    )


def rsid_column(df: narwhals.LazyFrame) -> narwhals.Expr:
    """The rsID of a study, or a column of nulls for a study that does not report one, so
    that studies remain vertically concatenable."""
    if GWASLAB_RSID_COL in df.collect_schema():
        return narwhals.col(GWASLAB_RSID_COL).cast(narwhals.String())
    return narwhals.lit(None, dtype=narwhals.String()).alias(GWASLAB_RSID_COL)


def prepare_source(
    df: narwhals.LazyFrame,
    sample_info: SampleInfo,
    keep_rsid: bool,
) -> narwhals.LazyFrame:
    """Reduce a study to its variant key and the quantities the meta analysis sums."""
    beta = narwhals.col(GWASLAB_BETA_COL).cast(narwhals.Float64())
    weight = 1 / narwhals.col(GWASLAB_SE_COL).cast(narwhals.Float64()) ** 2
    columns = [
        *[narwhals.col(c) for c in _VARIANT_KEY_COLS],
        weight.alias(_WEIGHT_COL),
        (beta * weight).alias(_WEIGHTED_BETA_COL),
        narwhals.lit(sample_info.effective_sample_size()).alias(
            _SOURCE_SAMPLE_SIZE_COL
        ),
    ]
    if keep_rsid:
        columns.append(rsid_column(df))
    return df.select(columns)


def combine_prepared_sources(
    dfs: Sequence[narwhals.LazyFrame], keep_rsid: bool
) -> narwhals.LazyFrame:
    """Combine prepared studies into the meta analysis by inverse variance weighting.

    Stacking the studies and grouping by variant, rather than joining them pairwise,
    keeps this to a single streaming pass whose memory is set by the number of distinct
    output variants rather than by the combined size of the studies.
    """
    aggregations = [
        narwhals.col(_WEIGHT_COL).sum().alias(_WEIGHT_COL),
        narwhals.col(_WEIGHTED_BETA_COL).sum().alias(_WEIGHTED_BETA_COL),
        narwhals.col(_SOURCE_SAMPLE_SIZE_COL).sum().alias(GWASLAB_SAMPLE_SIZE_COLUMN),
        narwhals.col(_WEIGHT_COL).count().alias(META_ANALYSIS_NUM_STUDIES_COL),
    ]
    if keep_rsid:
        # Any study reporting the variant gives the same rsID, so take whichever of them
        # is not null.
        aggregations.append(
            narwhals.col(GWASLAB_RSID_COL).min().alias(GWASLAB_RSID_COL)
        )
    aggregated = (
        narwhals.concat(dfs, how="vertical")
        .group_by(*_VARIANT_KEY_COLS)
        .agg(aggregations)
    )

    columns = [
        narwhals.col(GWASLAB_CHROM_COL),
        narwhals.col(GWASLAB_POS_COL),
        narwhals.col(GWASLAB_EFFECT_ALLELE_COL),
        narwhals.col(GWASLAB_NON_EFFECT_ALLELE_COL),
    ]
    if keep_rsid:
        columns.append(narwhals.col(GWASLAB_RSID_COL))
    columns.extend(
        [
            (narwhals.col(_WEIGHTED_BETA_COL) / narwhals.col(_WEIGHT_COL)).alias(
                GWASLAB_BETA_COL
            ),
            (1 / narwhals.col(_WEIGHT_COL)).sqrt().alias(GWASLAB_SE_COL),
            narwhals.col(GWASLAB_SAMPLE_SIZE_COLUMN).cast(narwhals.Int64()),
            narwhals.col(META_ANALYSIS_NUM_STUDIES_COL).cast(narwhals.Int64()),
        ]
    )
    # Grouping does not preserve any particular row order, and the build system's traces
    # hash the produced file, so sort on the full variant key to keep rebuilds of
    # unchanged inputs identical.
    return aggregated.select(columns).sort(_VARIANT_KEY_COLS)


def report_output_size(
    df: narwhals.LazyFrame,
) -> None:
    l = df.select(narwhals.len()).collect(engine="streaming").item()
    logger.debug(f"Final meta-analysis has {l} variants")


def report_study_count_distribution(df: narwhals.LazyFrame) -> None:
    """How many variants each number of contributing studies accounts for. Variants
    backed by a single study are the lowest powered part of a union meta analysis, so
    their share is worth seeing."""
    result = (
        df.group_by(META_ANALYSIS_NUM_STUDIES_COL)
        .agg(narwhals.len().alias("num_variants"))
        .sort(META_ANALYSIS_NUM_STUDIES_COL)
        .collect(engine="streaming")
    )
    logger.debug(f"Variants by number of contributing studies:\n{result}")
