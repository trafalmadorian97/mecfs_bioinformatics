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


def _acat_combine(p_values: pl.Series) -> float:
    """Combine k p-values via the Aggregated Cauchy Association Test (ACAT) with equal weights.

    Computes the equal-weight Cauchy test statistic
        T = (1 / k) * sum_{i=1..k} w_term(p_i)
    where k is the length of p_values, and w_term(p_i) is tan((0.5 - p_i) * pi) for
    p_i >= 1e-15 and its asymptotic form 1 / (pi * p_i) for 0 < p_i < 1e-15
    (Liu et al., 2019). The combined tail probability is evaluated as
        p_acat = atan2(1, T) / pi
    which equals 0.5 - atan(T) / pi without catastrophic cancellation for large positive T.

    Boundary values: returns 0.0 if any p_i is 0.0, and 1.0 if all p_i are 1.0.
    Individual p_i equal to 1.0 alongside non-unit p-values (common in discrete tests
    such as Fisher's exact test) are replaced by 1 - 1 / k so the pole at p = 1 does
    not cancel significant signals in other models.
    """
    k = len(p_values)
    assert k > 0, f"ACAT requires at least one p-value, got {k}"
    assert p_values.null_count() == 0, (
        f"ACAT p-values must not contain nulls, found {p_values.null_count()}"
    )

    vals: list[float] = p_values.cast(pl.Float64).to_list()
    for p in vals:
        assert not math.isnan(p) and 0.0 <= p <= 1.0, (
            f"ACAT p-values must be finite and in [0, 1], got {p}"
        )

    if any(p == 0.0 for p in vals):
        return 0.0
    if all(p == 1.0 for p in vals):
        return 1.0

    unit_replacement = 1.0 - 1.0 / k
    total = 0.0
    for p in vals:
        p_adj = unit_replacement if p == 1.0 else p
        if p_adj < _ACAT_SMALL_P_THRESHOLD:
            total += 1.0 / (math.pi * p_adj)
        else:
            total += math.tan((0.5 - p_adj) * math.pi)

    t = total / k
    return math.atan2(1.0, t) / math.pi


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

        result = (
            pl_df.group_by(self.group_by)
            .agg(
                pl.col(self.model_col).sort().str.join(", ").alias("models_used"),
                pl.col(self.model_col).count().alias("n_models"),
                pl.col(self.p_value_col)
                .map_batches(
                    _acat_combine,
                    return_dtype=pl.Float64,
                    returns_scalar=True,
                )
                .alias(ACAT_P_COL),
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
