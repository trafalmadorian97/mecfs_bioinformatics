"""Display table comparing two priors' annotation coefficients.

Each prior is explained by one or more annotation weights tables in the
RidgeAnnotationWeightsTask schema: e.g. a ridge surrogate of the precomputed
PolyFun prior, or the per-parity tau of a trait-specific L2-regularized S-LDSC
prior. A prior with several tables (one per chromosome parity) is summarized by
the per-annotation mean of its coefficients.

The coefficients compared are gamma_raw, the per-unit-annotation change in
per-SNP heritability. Both priors model per-SNP heritability linearly in the same
annotations, so for a given annotation their gamma_raw differ only by the
priors' overall scale, which carries no meaning here. Each prior's coefficients
are therefore divided by their largest absolute value, so the largest-magnitude
coefficient of each prior is +1 or -1.

gamma_standardized is not used: the two fits standardize on different scales
(per-SNP annotation SD for the ridge surrogate, LD-score SD for tau), so it is
not comparable across priors for a given annotation.
"""

from pathlib import Path, PurePath

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_RAW_COL,
    load_annotation_weights,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    float_column_names,
    write_parquet_table,
)
from mecfs_bio.build_system.wf.base_wf import WF

GAMMA_EXT_COL = "gamma_ext"
GAMMA_INT_COL = "gamma_int"


@frozen(slots=True)
class AnnotationWeightsComparisonTask(Task):
    """Write the external-vs-internal prior annotation coefficient table.

    Rows are sorted by the larger of the two normalized coefficient magnitudes,
    descending, so annotations either prior leans on come first.
    """

    meta: Meta
    external_weights_tasks: tuple[Task, ...]
    internal_weights_tasks: tuple[Task, ...]

    def __attrs_post_init__(self) -> None:
        assert self.external_weights_tasks, "Need at least one external weights table"
        assert self.internal_weights_tasks, "Need at least one internal weights table"

    @property
    def deps(self) -> list["Task"]:
        return [*self.external_weights_tasks, *self.internal_weights_tasks]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        external = _normalized_mean_gamma(
            [load_annotation_weights(fetch, t) for t in self.external_weights_tasks],
            GAMMA_EXT_COL,
        )
        internal = _normalized_mean_gamma(
            [load_annotation_weights(fetch, t) for t in self.internal_weights_tasks],
            GAMMA_INT_COL,
        )
        table = _comparison_table(external=external, internal=internal)
        out_path = scratch_dir / "annotation_weights_comparison.parquet"
        arrow_table = table.to_arrow()
        write_parquet_table(
            table=arrow_table,
            out_path=out_path,
            compression="zstd",
            compression_level=None,
            byte_stream_split_columns=float_column_names(arrow_table),
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        external_weights_tasks: tuple[Task, ...],
        internal_weights_tasks: tuple[Task, ...],
    ) -> "AnnotationWeightsComparisonTask":
        """Trait and project come from the first internal weights table, which
        belongs to the trait the internal prior was fit on."""
        source_meta = internal_weights_tasks[0].meta
        if not isinstance(source_meta, (ResultTableMeta, ResultDirectoryMeta)):
            raise ValueError(f"Unknown meta for internal weights task: {source_meta}")
        meta = ResultTableMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            external_weights_tasks=external_weights_tasks,
            internal_weights_tasks=internal_weights_tasks,
        )


def _normalized_mean_gamma(
    weights_tables: list[pl.DataFrame], gamma_col: str
) -> pl.DataFrame:
    """Per-annotation mean gamma_raw across the tables, divided by its largest
    absolute value. Keeps the family column."""
    annotations = [set(t[ANNOTATION_COL].to_list()) for t in weights_tables]
    assert all(a == annotations[0] for a in annotations), (
        "Weights tables of one prior cover different annotations"
    )
    mean = (
        pl.concat(
            [
                t.select(ANNOTATION_COL, FAMILY_COL, GAMMA_RAW_COL)
                for t in weights_tables
            ]
        )
        .group_by(ANNOTATION_COL, FAMILY_COL)
        .agg(pl.col(GAMMA_RAW_COL).mean().alias(gamma_col))
    )
    assert mean.height == len(annotations[0]), (
        "An annotation has different families across one prior's weights tables"
    )
    scale = float(np.abs(mean[gamma_col].to_numpy()).max())
    assert scale > 0.0, "Every coefficient of the prior is zero"
    return mean.with_columns(pl.col(gamma_col) / scale)


def _comparison_table(external: pl.DataFrame, internal: pl.DataFrame) -> pl.DataFrame:
    assert set(external[ANNOTATION_COL]) == set(internal[ANNOTATION_COL]), (
        "The two priors are explained by different annotations"
    )
    joined = external.join(internal, on=[ANNOTATION_COL, FAMILY_COL], how="inner")
    assert joined.height == external.height, (
        "The two priors assign some annotation to different families"
    )
    return joined.select(FAMILY_COL, ANNOTATION_COL, GAMMA_EXT_COL, GAMMA_INT_COL).sort(
        pl.max_horizontal(pl.col(GAMMA_EXT_COL).abs(), pl.col(GAMMA_INT_COL).abs()),
        descending=True,
    )
