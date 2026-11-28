"""
Task to combine p-values across collapsing models via the Aggregated Cauchy Association Test (ACAT).

Liu, Yaowu, et al. "ACAT: a fast and powerful p value combination method for rare-variant
analysis in sequencing studies." The American Journal of Human Genetics 104.3 (2019): 410-421.
"""

import math
from pathlib import Path, PurePath
from typing import Sequence

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

ACAT_P_COL = "acat_p"

_ACAT_SMALL_P_THRESHOLD = 1e-15


def acat_combine_expr(p: pl.Expr) -> pl.Expr:
    """Equal-weight ACAT over the k p-values in each group (Liu et al., 2019).

    T = mean_i c(p_i), where c(p) = tan((0.5 - p) * pi) for p >= 1e-15 and its
    asymptotic form 1 / (pi * p) below that. The combined p-value is atan2(1, T) / pi,
    which equals 0.5 - atan(T) / pi without cancellation for large T.
    Returns 0 if any p_i is 0 and 1 if all p_i are 1. A p_i of 1 alongside other
    p-values is replaced by 1 - 1/k so the pole at p = 1 does not cancel real signal.
    """
    p_adj = pl.when(p == 1.0).then(1.0 - 1.0 / p.len()).otherwise(p)
    cauchy_term = (
        pl.when(p_adj < _ACAT_SMALL_P_THRESHOLD)
        .then(1.0 / (math.pi * p_adj))
        .otherwise(((0.5 - p_adj) * math.pi).tan())
    )
    combined = pl.arctan2(pl.lit(1.0), cauchy_term.mean()) / math.pi
    return (
        pl.when((p == 0.0).any())
        .then(0.0)
        .when((p == 1.0).all())
        .then(1.0)
        .otherwise(combined)
    )


def _assert_valid_p_values(p_values: pl.Series) -> None:
    assert p_values.null_count() == 0, (
        f"ACAT p-values must not contain nulls, found {p_values.null_count()}"
    )
    invalid = p_values.filter(p_values.is_nan() | ~p_values.is_between(0.0, 1.0))
    assert invalid.is_empty(), (
        f"ACAT p-values must be in [0, 1], got {invalid.head(5).to_list()}"
    )


@frozen(slots=True)
class AcatTask(Task):
    source_task: Task
    meta: Meta
    group_by: list[str]
    p_value_col: str
    model_col: str
    excluded_models: list[str]

    @property
    def deps(self) -> list["Task"]:
        return [self.source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        asset = fetch(self.source_task.asset_id)
        lf = scan_dataframe_asset(
            asset, self.source_task.meta, parquet_backend="polars"
        )
        pl_df = pl.from_dataframe(lf.collect())

        if self.excluded_models:
            pl_df = pl_df.filter(~pl.col(self.model_col).is_in(self.excluded_models))

        _assert_valid_p_values(pl_df[self.p_value_col].cast(pl.Float64))
        result = (
            pl_df.group_by(self.group_by)
            .agg(
                pl.col(self.model_col).sort().str.join(", ").alias("models_used"),
                pl.col(self.model_col).count().alias("n_models"),
                acat_combine_expr(pl.col(self.p_value_col).cast(pl.Float64)).alias(
                    ACAT_P_COL
                ),
            )
            .sort(self.group_by)
        )

        out_path = scratch_dir / "acat_result.parquet"
        result.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        source_task: Task,
        asset_id: str,
        group_by: Sequence[str],
        p_value_col: str,
        model_col: str,
        excluded_models: Sequence[str] = (),
    ) -> "AcatTask":
        source_meta = source_task.meta
        assert isinstance(source_meta, (GWASSummaryDataFileMeta, FilteredGWASDataMeta))
        meta = FilteredGWASDataMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("processed"),
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(
            source_task=source_task,
            meta=meta,
            group_by=list(group_by),
            p_value_col=p_value_col,
            model_col=model_col,
            excluded_models=list(excluded_models),
        )
