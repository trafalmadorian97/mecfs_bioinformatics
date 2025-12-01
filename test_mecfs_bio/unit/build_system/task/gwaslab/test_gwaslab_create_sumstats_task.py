import pickle
from pathlib import Path

import gwaslab

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_gwaslab_sumstats(
    tmp_path: Path,
):
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = GWASLabCreateSumstatsTask(
        asset_id=AssetId("sumstats_task"),
        basic_check=True,
        df_source_task=FakeTask(
            meta=FilteredGWASDataMeta(
                AssetId("input"),
                "dummy_trait",
                "dummy_project",
                "dummy_Dir",
                read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=" ")),
            )
        ),
        genome_build="38",
    )

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(
            Path("test_mecfs_bio/unit/build_system/task/dummy_data.regenie")
        )

    asset_result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    with open(asset_result.path, "rb") as f:
        loaded = pickle.load(
            f,
        )
    assert isinstance(loaded, gwaslab.Sumstats)
    scratch_loc_2 = tmp_path / "scratch_2"
    scratch_loc_2.mkdir(exist_ok=True, parents=True)

    def fetch_2(asset_id: AssetId) -> Asset:
        return asset_result

    task_2 = GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=task, asset_id="table_task", sub_dir="processed"
    )
    asset_2 = task_2.execute(scratch_dir=scratch_loc_2, wf=task_2, fetch=fetch_2)
    scan_dataframe_asset(asset=asset_2, meta=task_2.meta)
