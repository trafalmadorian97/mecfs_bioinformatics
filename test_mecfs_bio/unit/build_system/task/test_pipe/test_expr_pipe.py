import narwhals
import numpy as np
import polars as pl

from mecfs_bio.build_system.task.pipes.expr_pipe import ExprPipe


def test_expr_pipe():
    df = pl.DataFrame({"col": [1, 4, 8]})
    result = (
        ExprPipe((2 * narwhals.col("col")).alias("col"))
        .process(narwhals.from_native(df).lazy())
        .collect()
        .to_polars()
    )
    np.testing.assert_allclose(result["col"].to_numpy(), np.array([2, 8, 16]))
