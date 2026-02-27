from pathlib import Path

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.combine_gene_lists_task import (
    CombineGeneListsTask,
    SrcGeneList,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_combine_gene_lists(tmp_path: Path):
    df_1_loc = tmp_path / "df_1.csv"
    df_2_loc = tmp_path / "df_2.csv"
    df_1 = pd.DataFrame({"gene1": ["A", "B", "C"]})
    df_1.to_csv(df_1_loc, index=False)
    df_2 = pd.DataFrame({"gene2": ["B", "C", "D"]})
    df_2.to_csv(df_2_loc, index=False)

    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)

    task = CombineGeneListsTask(
        meta=SimpleFileMeta(AssetId("combine")),
        src_gene_lists=[
            SrcGeneList(
                task=FakeTask(
                    SimpleFileMeta(
                        AssetId("df1"), DataFrameReadSpec(DataFrameTextFormat(","))
                    ),
                ),
                name="df1",
                ensemble_id_column="gene1",
            ),
            SrcGeneList(
                task=FakeTask(
                    SimpleFileMeta(
                        AssetId("df2"), DataFrameReadSpec(DataFrameTextFormat(","))
                    ),
                ),
                name="df2",
                ensemble_id_column="gene2",
            ),
        ],
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "df1":
            return FileAsset(df_1_loc)
        if asset_id == "df2":
            return FileAsset(df_2_loc)
        raise ValueError(f"unknown asset id {asset_id}")

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_csv(result.path)
    assert len(result_df) == 4
