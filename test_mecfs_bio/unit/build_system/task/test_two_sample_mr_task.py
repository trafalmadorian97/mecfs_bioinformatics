import pandas as pd
import pytest

from mecfs_bio.build_system.task.two_sample_mr_task import (
    TSM_BETA_COL,
    TSM_EFFECT_ALLELE_COL,
    TSM_OTHER_ALLELE_COL,
    TSM_OUTPUT_B_COL,
    TSM_RSID_COL,
    TSM_SE_COL,
    TwoSampleMRConfig,
    run_two_sample_mr,
)


def test_two_sample_mr_task():
    """
    Test that we get the correct result in dummy example of performing MR
    """
    dummy_exposure = pd.DataFrame(
        {
            TSM_RSID_COL: ["rs12345", "rs456"],
            TSM_BETA_COL: [0.5, 1],
            TSM_SE_COL: [0.01, 0.02],
            TSM_EFFECT_ALLELE_COL: ["A", "G"],
            TSM_OTHER_ALLELE_COL: ["C", "T"],
        }
    )

    dummy_outcome = pd.DataFrame(
        {
            TSM_RSID_COL: ["rs12345", "rs456"],
            TSM_BETA_COL: [1, 2.1],
            TSM_SE_COL: [0.01, 0.02],
            TSM_EFFECT_ALLELE_COL: ["A", "G"],
            TSM_OTHER_ALLELE_COL: ["C", "T"],
        }
    )
    tsm_result = run_two_sample_mr(
        exposure_df=dummy_exposure,
        outcome_df=dummy_outcome,
        config=TwoSampleMRConfig(clump_exposure_data=False),
    )
    assert tsm_result.result[TSM_OUTPUT_B_COL].item() == pytest.approx(2, abs=0.1)
