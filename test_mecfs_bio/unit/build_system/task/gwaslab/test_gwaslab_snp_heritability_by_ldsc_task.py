import numpy as np
import pandas as pd

from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    _drop_variants_with_degenerate_z,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_P_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ldsc_constants import LDSC_Z_COL

_ID_COL = "SNP"


def test_drops_zero_se_variants():
    # SE == 0 with BETA == 0 (odds ratio rounded to 1.00) -> Z would be NaN;
    # SE == 0 with BETA != 0 (underflowed p-value) -> Z would be inf.
    data = pd.DataFrame(
        {
            _ID_COL: ["a", "b", "c", "d"],
            GWASLAB_BETA_COL: [0.1, 0.0, 0.2, -0.3],
            GWASLAB_SE_COL: [0.05, 0.0, 0.0, 0.04],
        }
    )
    result = _drop_variants_with_degenerate_z(data)
    assert result[_ID_COL].tolist() == ["a", "d"]


def test_keeps_all_finite_se_variants():
    data = pd.DataFrame(
        {
            _ID_COL: ["a", "b"],
            GWASLAB_BETA_COL: [0.1, -0.2],
            GWASLAB_SE_COL: [0.05, 0.03],
        }
    )
    result = _drop_variants_with_degenerate_z(data)
    assert result[_ID_COL].tolist() == ["a", "b"]


def test_prefers_existing_z_column():
    data = pd.DataFrame(
        {
            _ID_COL: ["a", "b", "c"],
            LDSC_Z_COL: [1.5, np.inf, np.nan],
            # BETA/SE would keep every row; the Z column must take precedence.
            GWASLAB_BETA_COL: [0.1, 0.2, 0.3],
            GWASLAB_SE_COL: [0.05, 0.05, 0.05],
        }
    )
    result = _drop_variants_with_degenerate_z(data)
    assert result[_ID_COL].tolist() == ["a"]


def test_noop_without_z_or_beta_se_columns():
    data = pd.DataFrame({_ID_COL: ["a", "b"], GWASLAB_P_COL: [0.1, 0.2]})
    result = _drop_variants_with_degenerate_z(data)
    assert result[_ID_COL].tolist() == ["a", "b"]
