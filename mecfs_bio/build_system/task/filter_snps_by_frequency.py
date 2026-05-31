"""
Task to filter SNPs in a GWAS based on a minimum minor allele frequency.
"""

from pathlib import Path, PurePath

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

import narwhals as nw


@frozen
class FilterSNPsFrequencyTask(Task):
    """
    Task to only keep variants where the minor allele frequency
    is greater than or equal to the given freq_thresh.
    """

    raw_gwas_task: Task
    meta: FilteredGWASDataMeta | SimpleFileMeta
    allele_freq_col: str = "effect_allele_frequency"
    freq_thresh: float = 0.05

    @property
    def deps(self) -> list["Task"]:
        return [self.raw_gwas_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        logger.debug(f"Filtering to common variants (MAF >= {self.freq_thresh})")
        df = scan_dataframe_asset(
            asset=fetch(self.raw_gwas_task.asset_id),
            meta=self.raw_gwas_task.meta,
        )
        result = df.filter(
            nw.min_horizontal(
                nw.col(self.allele_freq_col), 1 - nw.col(self.allele_freq_col)
            )
            >= self.freq_thresh
        )
        target_path = scratch_dir / "tmp.parqet"
        result.sink_parquet(target_path)
        return FileAsset(target_path)

    @classmethod
    def create(
        cls,
        id: str,
        raw_gwas_task: Task,
        allele_freq_col: str = "effect_allele_frequency",
        freq_thresh: float = 0.05,
    ) -> Task:
        source_meta = raw_gwas_task.meta
        assert isinstance(source_meta, GWASSummaryDataFileMeta)
        meta = FilteredGWASDataMeta(
            id=AssetId(id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("processed"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(format=DataFrameParquetFormat()),
        )
        return cls(
            raw_gwas_task=raw_gwas_task,
            meta=meta,
            allele_freq_col=allele_freq_col,
            freq_thresh=freq_thresh,
        )
