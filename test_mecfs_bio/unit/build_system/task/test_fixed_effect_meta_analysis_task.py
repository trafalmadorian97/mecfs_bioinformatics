from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fixed_effect_meta_analysis_task import (
    CaseControlSampleInfo,
    FixedEffectsMetaAnalysisTask,
    GwasSource,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)


def test_fixed_effect_meta_analysis_task(tmp_path: Path):
    df_1 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [1],
            GWASLAB_SE_COL: [0.5],
        }
    )

    df_2 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [1],
            GWASLAB_SE_COL: [0.25],
        }
    )
    df_3 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [10],
            GWASLAB_SE_COL: [10],
        }
    )
    df_1_path = tmp_path / "df_1.parquet"
    df_1.to_parquet(df_1_path)
    df_2_path = tmp_path / "df_2.parquet"
    df_2.to_parquet(df_2_path)
    df_3_path = tmp_path / "df_3.parquet"
    df_3.to_parquet(df_3_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    task = FixedEffectsMetaAnalysisTask(
        meta=SimpleFileMeta("meta_analysis"),
        sources=[
            GwasSource(
                task=FakeTask(
                    SimpleFileMeta(
                        "df1", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                sample_info=CaseControlSampleInfo(controls=10, cases=10),
            ),
            GwasSource(
                task=FakeTask(
                    SimpleFileMeta(
                        "df2", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                sample_info=CaseControlSampleInfo(controls=10, cases=10),
            ),
            GwasSource(
                task=FakeTask(
                    SimpleFileMeta(
                        "df3", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                sample_info=CaseControlSampleInfo(controls=10, cases=10),
            ),
        ],
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "df1":
            return FileAsset(df_1_path)
        if asset_id == "df2":
            return FileAsset(df_2_path)
        if asset_id == "df3":
            return FileAsset(df_3_path)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    df = pd.read_parquet(result.path)
    assert df["SE"].item() == pytest.approx(0.2235, abs=0.01)
