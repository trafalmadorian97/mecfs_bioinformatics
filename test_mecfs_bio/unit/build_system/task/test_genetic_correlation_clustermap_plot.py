from pathlib import Path

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.data_manipulation.xr_data.pipes.xr_identity import (
    XRIdentityPipe,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genetic_correlation_clustermap_task import (
    GeneticCorrelationClustermapTask,
    GeneticCorrSource,
    RGWithAsterix,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def test_genetic_correlation_clustermap(tmp_path: Path):
    df_loc = tmp_path / "df.csv"
    df = pd.DataFrame(
        {
            "p1": ["multisite_pain", "multisite_pain", "schizophrenia"],
            "p2": ["schizophrenia", "asthma", "asthma"],
            "rg": [0.0780454039363654, 0.2669804689854996, 0.0582518935695738],
            "p": [0.0005428096193185, 2.0149832056929464e-25, 0.0107988881808604],
        }
    )
    df.to_csv(df_loc, index=False)
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = GeneticCorrelationClustermapTask(
        meta=SimpleFileMeta(AssetId("plot")),
        genetic_corr_source=GeneticCorrSource(
            task=FakeTask(
                SimpleFileMeta(
                    AssetId(
                        "corr_df",
                    ),
                    read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
                )
            )
        ),
        plot_options=RGWithAsterix(),
        xr_pipe=XRIdentityPipe(),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "corr_df":
            return FileAsset(df_loc)
        raise NotImplementedError()

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
