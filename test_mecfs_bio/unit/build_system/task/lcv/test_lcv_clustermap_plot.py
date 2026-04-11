from pathlib import Path

import pandas as pd

from mecfs_bio.asset_generator.lcv_asset_generator import (
    DOWNSTREAM_TRAIT_COL,
    UPSTREAM_TRAIT_COL,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.lcv.lcv_clustermap import (
    GCPWithAsterisk,
    LCVClustermapTask,
    LCVSource,
)
from mecfs_bio.build_system.task.lcv.lcv_core import (
    LCV_MEAN_GCP_COL,
    LCV_PVAL_ZERO_COL,
    LCV_RHO_EST_COL,
)
from mecfs_bio.build_system.task.xr_pipes.xr_identity import XRIdentityPipe
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_lcv_clustermap(tmp_path: Path):
    df_loc = tmp_path / "lcv_agg.parquet"
    df = pd.DataFrame(
        {
            UPSTREAM_TRAIT_COL: ["mecfs", "mecfs", "crp", "crp"],
            DOWNSTREAM_TRAIT_COL: [
                "crp",
                "schizophrenia",
                "mecfs",
                "schizophrenia",
            ],
            LCV_RHO_EST_COL: [0.12, 0.05, 0.18, -0.03],
            LCV_PVAL_ZERO_COL: [1e-10, 0.4, 1e-3, 0.8],
            LCV_MEAN_GCP_COL: [0.72, 0.05, -0.30, 0.10],
        }
    )
    df.to_parquet(df_loc, index=False)

    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)

    task = LCVClustermapTask(
        meta=SimpleFileMeta(AssetId("lcv_plot")),
        xr_pipe=XRIdentityPipe(),
        source=LCVSource(
            task=FakeTask(
                SimpleFileMeta(
                    AssetId("lcv_df"),
                    read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                )
            ),
        ),
        plot_options=GCPWithAsterisk(),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "lcv_df":
            return FileAsset(df_loc)
        raise NotImplementedError()

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    assert result.path.exists()
    assert result.path.stat().st_size > 0
