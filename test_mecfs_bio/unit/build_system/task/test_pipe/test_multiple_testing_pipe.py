import narwhals
import numpy as np
import polars as pl

from mecfs_bio.build_system.task.multiple_testing_table_task import REJECT_NULL_LABEL
from mecfs_bio.build_system.task.pipes.multiple_testing_pipe import MultipleTestingPipe


def test_multiple_testing_pipe():
    df = pl.DataFrame({"p": [0.000001, 0.05, 1]})
    result = (
        MultipleTestingPipe("p")
        .process(narwhals.from_native(df).lazy())
        .collect()
        .to_polars()
    )
    np.testing.assert_allclose(
        result[REJECT_NULL_LABEL].to_numpy(), np.array([True, False, False])
    )
