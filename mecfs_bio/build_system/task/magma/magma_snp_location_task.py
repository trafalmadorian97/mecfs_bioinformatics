from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.assign_rsids_via_snp151_task import create_new_meta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_P_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.build_system.task.pipes.compute_p_pipe import (
    ComputePIfNeededPipe,
    ComputePPipe,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class MagmaSNPFileTask(Task):
    gwas_parquet_with_rsid_task: Task
    extra_columns_to_output: list[str]
    _meta: Meta
    pipes: Sequence[DataProcessingPipe]

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def source_meta(self) -> Meta:
        return self.gwas_parquet_with_rsid_task.meta

    @property
    def source_id(self) -> AssetId:
        return self.source_meta.asset_id

    @property
    def deps(self) -> list["Task"]:
        return [self.gwas_parquet_with_rsid_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.source_id)
        gwas_data = scan_dataframe_asset(source_asset, meta=self.source_meta)
        for pipe in self.pipes:
            gwas_data = pipe.process(gwas_data)
        gwas_data = gwas_data.unique(
            subset=GWASLAB_RSID_COL, order_by=[GWASLAB_CHROM_COL, GWASLAB_POS_COL]
        )
        out_path = scratch_dir / "output"
        gwas_data.sort([GWASLAB_CHROM_COL, GWASLAB_POS_COL]).select(
            GWASLAB_RSID_COL, *self.extra_columns_to_output
        ).collect().to_polars().write_csv(out_path, include_header=False, separator=" ")
        return FileAsset(
            out_path,
        )

    @classmethod
    def create_for_magma_snp_pos_file(
        cls,
        gwas_parquet_with_rsids_task: Task,
        asset_id: str,
        pipes: list[DataProcessingPipe] | None = None,
    ):
        if pipes is None:
            pipes = []
        extra_cols = [GWASLAB_CHROM_COL, GWASLAB_POS_COL]
        source_meta = gwas_parquet_with_rsids_task.meta
        meta = create_new_meta(
            source_meta,
            asset_id=asset_id,
            format=DataFrameTextFormat(
                separator=" ",
                has_header=False,
                column_names=[GWASLAB_RSID_COL] + extra_cols,
            ),
            extension=".id.chr.genpos.txt",
        )
        return cls(
            meta=meta,
            gwas_parquet_with_rsid_task=gwas_parquet_with_rsids_task,
            extra_columns_to_output=extra_cols,
            pipes=pipes,
        )

    @classmethod
    def create_for_magma_snp_p_value_file(
        cls,
        gwas_parquet_with_rsids_task: Task,
        asset_id: str,
        pipes: list[DataProcessingPipe] | None = None,
    ):
        if pipes is None:
            pipes = []
        pipes.append(ComputePPipe())
        extra_cols = [GWASLAB_P_COL]
        source_meta = gwas_parquet_with_rsids_task.meta
        meta = create_new_meta(
            source_meta,
            asset_id=asset_id,
            format=DataFrameTextFormat(
                separator=" ",
                has_header=False,
                column_names=[GWASLAB_RSID_COL] + extra_cols,
            ),
            extension=".id.p.txt",
        )
        return cls(
            meta=meta,
            gwas_parquet_with_rsid_task=gwas_parquet_with_rsids_task,
            extra_columns_to_output=extra_cols,
            pipes=pipes,
        )

    @classmethod
    def create_for_magma_snp_p_value_file_precomputed_p(
        cls,
        gwas_parquet_with_rsids_task: Task,
        asset_id: str,
    ):
        """
        As above, but assume the data already contains p values.
        """
        extra_cols = [GWASLAB_P_COL]
        source_meta = gwas_parquet_with_rsids_task.meta
        meta = create_new_meta(
            source_meta,
            asset_id=asset_id,
            format=DataFrameTextFormat(
                separator=" ",
                has_header=False,
                column_names=[GWASLAB_RSID_COL] + extra_cols,
            ),
            extension=".id.p.txt",
        )
        return cls(
            meta=meta,
            gwas_parquet_with_rsid_task=gwas_parquet_with_rsids_task,
            extra_columns_to_output=extra_cols,
            pipes=[],
        )

    @classmethod
    def create_for_magma_snp_p_value_file_compute_if_needed(
        cls,
        gwas_parquet_with_rsids_task: Task,
        asset_id: str,
        pipes: list[DataProcessingPipe] | None = None,
    ):
        if pipes is None:
            pipes = []
        pipes.append(ComputePIfNeededPipe())
        extra_cols = [GWASLAB_P_COL]
        source_meta = gwas_parquet_with_rsids_task.meta
        meta = create_new_meta(
            source_meta,
            asset_id=asset_id,
            format=DataFrameTextFormat(
                separator=" ",
                has_header=False,
                column_names=[GWASLAB_RSID_COL] + extra_cols,
            ),
            extension=".id.p.txt",
        )
        return cls(
            meta=meta,
            gwas_parquet_with_rsid_task=gwas_parquet_with_rsids_task,
            extra_columns_to_output=extra_cols,
            pipes=pipes,
        )
