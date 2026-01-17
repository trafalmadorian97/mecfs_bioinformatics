import narwhals
import numpy as np
import pandas as pd

from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_P_COL,
    GWASLAB_SE_COL,
)


def test_compute_se_pipe():
    dummy_data = pd.DataFrame({GWASLAB_BETA_COL: [0.0388], GWASLAB_P_COL: [0.45056]})
    nw_data = narwhals.from_native(dummy_data).lazy()
    pipe = ComputeSEPipe()
    result = pipe.process(nw_data).collect().to_pandas()
    np.testing.assert_allclose(
        result[GWASLAB_SE_COL].item(), 0.05142602, rtol=0.0001, atol=0.0001
    )


def test_compute_se_pipe_handles_small_pvalues():
    """
    Verify that our pipe to compute standard errors is able to handle small p values
    """

    dummy_data = pd.DataFrame(
        {GWASLAB_BETA_COL: [0.0388], GWASLAB_P_COL: [6.500000e-18]}
    )
    nw_data = narwhals.from_native(dummy_data).lazy()
    pipe = ComputeSEPipe()
    result = pipe.process(nw_data).collect().to_pandas()
    assert result[GWASLAB_SE_COL].item() > 0.0
