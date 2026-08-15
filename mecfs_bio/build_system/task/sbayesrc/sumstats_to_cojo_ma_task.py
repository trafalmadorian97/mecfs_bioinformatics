"""
Task converting a GWAS summary statistics asset into GCTB/COJO .ma format.

A .ma file is GCTB's plain-text input format for summary statistics.
"""

from pathlib import Path, PurePath

import narwhals as nw
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import COJO_MA_COLUMNS
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

_SNP_COL, _A1_COL, _A2_COL, _FREQ_COL, _B_COL, _SE_COL, _P_COL, _N_COL = COJO_MA_COLUMNS


@frozen
class SumstatsToCojoMaTask(Task):
    """
    Reads a tabular summary statistics file and writes it in GCTB/COJO .ma format.
    """

    sumstats_task: Task
    meta: FilteredGWASDataMeta

    @property
    def deps(self) -> list["Task"]:
        return [self.sumstats_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        source_lf = scan_dataframe_asset(
            asset=fetch(self.sumstats_task.asset_id), meta=self.sumstats_task.meta
        )
        selection = source_lf.select(
            nw.col(GWASLAB_RSID_COL).alias(_SNP_COL),
            nw.col(GWASLAB_EFFECT_ALLELE_COL).alias(_A1_COL),
            nw.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(_A2_COL),
            nw.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).alias(_FREQ_COL),
            nw.col(GWASLAB_BETA_COL).alias(_B_COL),
            nw.col(GWASLAB_SE_COL).alias(_SE_COL),
            nw.col(GWASLAB_P_COL).alias(_P_COL),
            nw.col(GWASLAB_SAMPLE_SIZE_COLUMN).alias(_N_COL),
        )
        result_df = selection.collect().to_polars()
        target_path = scratch_dir / f"{self.meta.id}.ma"
        result_df.write_csv(target_path, separator="\t")
        return FileAsset(target_path)

    @classmethod
    def create(cls, id: str, sumstats_task: Task) -> "SumstatsToCojoMaTask":
        source_meta = sumstats_task.meta
        assert isinstance(source_meta, (GWASSummaryDataFileMeta, FilteredGWASDataMeta))
        assert isinstance(source_meta.read_spec, DataFrameReadSpec)
        meta = FilteredGWASDataMeta(
            id=AssetId(id),
            trait=source_meta.trait,
            project=source_meta.project,
            sub_dir=PurePath("processed") / "gwfm" / "cojo_ma",
            extension=".ma",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
        )
        return cls(sumstats_task=sumstats_task, meta=meta)
