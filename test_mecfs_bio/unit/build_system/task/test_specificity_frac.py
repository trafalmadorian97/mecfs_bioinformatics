from pathlib import Path
from typing import Any

import pandas as pd

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.specificity_frac_task import (
    NORMALIZED_MEAN,
    PrepareSpecificityFraction,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

_dummy_df = pd.DataFrame(  # According to the fractional metric, Y is specific for A, while X is specifific for V
    {
        "cell_type": ["A", "A", "A", "A", "B", "B", "B", "B"],
        "cell": ["A1", "A1", "A2", "A2", "B1", "B1", "B2", "B2"],
        "gene": ["X", "Y", "X", "Y", "X", "Y", "X", "Y"],
        "count": [10, 1, 10, 50, 1, 10, 50, 10],
    }
)


def test_prepare_specificity_frac(tmp_path: Path):
    """
    Test computing the specificity of genes for cell types using the fractional approach
    """
    df_loc = tmp_path / "df_loc.parquet"
    _dummy_df.to_parquet(df_loc)
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    tsk = PrepareSpecificityFraction(
        meta=SimpleFileMeta("simple_fle"),
        long_count_df_task=FakeTask(
            SimpleFileMeta(
                "count_df", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            )
        ),
        cell_type_col="cell_type",
        count_col="count",
        gene_col="gene",
        cell_col="cell",
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "count_df":
            return FileAsset(df_loc)
        raise ValueError("unknown asset id")

    result = tsk.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    result_df = pd.read_parquet(result.path)

    pivoted: Any = result_df.pivot(
        index="cell_type", columns="gene", values=NORMALIZED_MEAN
    )
    assert float(pivoted.loc["A", "X"]) < float(pivoted.loc["B", "X"])
    assert float(pivoted.loc["A", "Y"]) > float(pivoted.loc["B", "Y"])
