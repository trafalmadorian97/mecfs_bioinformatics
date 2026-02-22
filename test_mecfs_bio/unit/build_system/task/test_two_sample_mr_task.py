from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.two_sample_mr_task import (
    MAIN_RESULT_DF_PATH,
    TSM_BETA_COL,
    TSM_EFFECT_ALLELE_COL,
    TSM_OTHER_ALLELE_COL,
    TSM_OUTPUT_B_COL,
    TSM_OUTPUT_EXPOSURE_COL,
    TSM_PHENOTYPE,
    TSM_RSID_COL,
    TSM_SE_COL,
    TwoSampleMRConfig,
    TwoSampleMRTask,
    run_two_sample_mr,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF

dummy_exposure = pd.DataFrame(
    {
        TSM_RSID_COL: ["rs12345", "rs456"],
        TSM_BETA_COL: [0.5, 1],
        TSM_SE_COL: [0.01, 0.02],
        TSM_EFFECT_ALLELE_COL: ["A", "G"],
        TSM_OTHER_ALLELE_COL: ["C", "T"],
    }
)

dummy_exposure_A = pd.DataFrame(
    {
        TSM_PHENOTYPE: ["A"] * 2,
        TSM_RSID_COL: ["rs12345", "rs456"],
        TSM_BETA_COL: [0.5, 1],
        TSM_SE_COL: [0.01, 0.02],
        TSM_EFFECT_ALLELE_COL: ["A", "G"],
        TSM_OTHER_ALLELE_COL: ["C", "T"],
    }
)

dummy_exposure_B = pd.DataFrame(
    {
        TSM_PHENOTYPE: ["B"] * 2,
        TSM_RSID_COL: ["rs789", "rs101112"],
        TSM_BETA_COL: [1, 2],
        TSM_SE_COL: [0.01, 0.02],
        TSM_EFFECT_ALLELE_COL: ["A", "G"],
        TSM_OTHER_ALLELE_COL: ["C", "T"],
    }
)

dummy_outcome_A = pd.DataFrame(
    {
        TSM_RSID_COL: ["rs12345", "rs456"],
        TSM_BETA_COL: [1, 2.1],
        TSM_SE_COL: [0.01, 0.02],
        TSM_EFFECT_ALLELE_COL: ["A", "G"],
        TSM_OTHER_ALLELE_COL: ["C", "T"],
    }
)

dummy_outcome_B = pd.DataFrame(
    {
        TSM_RSID_COL: ["rs789", "rs101112"],
        TSM_BETA_COL: [1, 2.1],
        TSM_SE_COL: [0.01, 0.02],
        TSM_EFFECT_ALLELE_COL: ["A", "G"],
        TSM_OTHER_ALLELE_COL: ["C", "T"],
    }
)


def test_two_sample_mr_func():
    """
    Test that we get the correct result in dummy example of performing MR
    """
    tsm_result = run_two_sample_mr(
        exposure_df=dummy_exposure,
        outcome_df=dummy_outcome_A,
        config=TwoSampleMRConfig(clump_exposure_data=None),
    )
    assert tsm_result.result[TSM_OUTPUT_B_COL].item() == pytest.approx(2, abs=0.1)


def test_two_sample_mr_with_multiple_exposure():
    """
    When we run MR with two input exposures, there should be a separate output for each exposure value
    """
    combined_exposure = pd.concat([dummy_exposure_A, dummy_exposure_B], axis=0)
    combined_outcome = pd.concat([dummy_outcome_A, dummy_outcome_B], axis=0)
    tsm_result = run_two_sample_mr(
        exposure_df=combined_exposure,
        outcome_df=combined_outcome,
        config=TwoSampleMRConfig(clump_exposure_data=None),
    )
    assert len(tsm_result.result[TSM_OUTPUT_EXPOSURE_COL].unique()) == 2


def test_two_sample_mr_task(tmp_path: Path):
    dummy_exposure_path = tmp_path / "dummy_exposure.csv"
    dummy_outcome_path = tmp_path / "dummy_outcome.csv"
    dummy_exposure.to_csv(dummy_exposure_path, index=False)
    dummy_outcome_A.to_csv(dummy_outcome_path, index=False)
    scratch_loc = tmp_path / "scratch"
    scratch_loc.mkdir(exist_ok=True, parents=True)
    task = TwoSampleMRTask(
        SimpleFileMeta(AssetId("mr_result")),
        outcome_data_task=FakeTask(
            SimpleFileMeta(
                AssetId("outcome_data"),
                read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
            ),
        ),
        exposure_data_task=FakeTask(
            SimpleFileMeta(
                AssetId("exposure_data"),
                read_spec=DataFrameReadSpec(DataFrameTextFormat(",")),
            ),
        ),
        config=TwoSampleMRConfig(clump_exposure_data=None),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "outcome_data":
            return FileAsset(dummy_outcome_path)
        if asset_id == "exposure_data":
            return FileAsset(dummy_exposure_path)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch_loc, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, DirectoryAsset)
    loaded_result = pd.read_csv((result.path) / MAIN_RESULT_DF_PATH)
    assert loaded_result[TSM_OUTPUT_B_COL].item() == pytest.approx(2, abs=0.1)
