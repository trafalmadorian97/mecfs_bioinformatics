"""Display table comparing an external and an internal prior's SUSIE runs
against the uniform run at one locus.

One row per variant in any of the three runs' credible sets, with the columns of
load_prior_comparison_variants: each run's credible-set number and PIP, and each
prior's lift. Sorted by the larger of the two prior PIPs, descending.
"""

from pathlib import Path, PurePath

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
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    float_column_names,
    write_parquet_table,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_variants import (
    CS_COLS,
    PriorComparisonRuns,
    load_prior_comparison_variants,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.polyfun_explain_display_columns import (
    PIP_PF_EXT_COL,
    PIP_PF_INT_COL,
)


@frozen(slots=True)
class PolyfunPriorComparisonTableTask(Task):
    """Write the external-vs-internal prior comparison display table."""

    meta: Meta
    runs: PriorComparisonRuns
    secondary_position: SecondaryPositionFromSnpid | None = None

    @property
    def deps(self) -> list["Task"]:
        return self.runs.tasks

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        variants = load_prior_comparison_variants(
            fetch, self.runs, self.secondary_position
        )
        table = variants.filter(
            pl.any_horizontal(pl.col(c).is_not_null() for c in CS_COLS)
        ).sort(pl.max_horizontal(PIP_PF_EXT_COL, PIP_PF_INT_COL), descending=True)
        out_path = scratch_dir / "prior_comparison_table.parquet"
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
        runs: PriorComparisonRuns,
        secondary_position: SecondaryPositionFromSnpid | None = None,
    ) -> "PolyfunPriorComparisonTableTask":
        source_meta = runs.susie_uniform_task.meta
        if not isinstance(source_meta, ResultDirectoryMeta):
            raise ValueError(f"Unknown meta for uniform susie task: {source_meta}")
        meta = ResultTableMeta(
            id=AssetId(asset_id),
            trait=source_meta.trait,
            project=source_meta.project,
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            sub_dir=PurePath("analysis"),
        )
        return cls(meta=meta, runs=runs, secondary_position=secondary_position)
