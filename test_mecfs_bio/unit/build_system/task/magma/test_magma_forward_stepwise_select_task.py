import numpy as np
import pandas as pd
import pytest

from mecfs_bio.build_system.task.magma.magma_forward_stepwise_select_task import (
    generate_mappers_from_wide_dataframe,
    generate_wide_dataframe,
)


def test_proportional_significance():
    """
    Test that a proportional significance value is consistent with manual calculation
    """
    df_marg = pd.read_csv(
        "test_mecfs_bio/unit/build_system/task/magma/dummy_marg_magma_output.txt",
        sep=r"\s+",
    )
    df_cond = pd.read_csv(
        "test_mecfs_bio/unit/build_system/task/magma/dummy_cond_magma_output.txt",
        sep=r"\s+",
    )
    df_wide = generate_wide_dataframe(df_cond, df_marg=df_marg)
    marg_dict, prop_sig_dict = generate_mappers_from_wide_dataframe(df_wide)
    assert prop_sig_dict[("Cluster242", "Cluster239")] == pytest.approx(
        np.log10(0.00048536) / np.log10(1.4544e-15)
    )
