import pandas as pd
import pytest

from mecfs_bio.build_system.task.two_sample_mr_task import (
    TSM_BETA_COL,
    TSM_EFFECT_ALLELE_COL,
    TSM_OTHER_ALLELE_COL,
    TSM_OUTPUT_B_COL,
    TSM_OUTPUT_EXPOSURE_COL,
    TSM_PHENOTYPE,
    TSM_RSID_COL,
    TSM_SE_COL,
    TwoSampleMRConfig,
    run_two_sample_mr,
)

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


def test_two_sample_mr_task():
    """
    Test that we get the correct result in dummy example of performing MR
    """
    tsm_result = run_two_sample_mr(
        exposure_df=dummy_exposure,
        outcome_df=dummy_outcome_A,
        config=TwoSampleMRConfig(clump_exposure_data=False),
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
        config=TwoSampleMRConfig(clump_exposure_data=False),
    )
    assert len(tsm_result.result[TSM_OUTPUT_EXPOSURE_COL].unique()) == 2
