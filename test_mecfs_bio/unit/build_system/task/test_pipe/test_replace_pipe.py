import narwhals
import numpy as np
import polars as pl

from mecfs_bio.build_system.task.pipes.replace_pipe import ReplaceStrictPipe


def test_replace_pipe():
    df = pl.DataFrame({"col": [float("nan"), 1.5]})
    result = (
        ReplaceStrictPipe(
            target_column="col",
            new_column_name="col",
            replace_mapping={float("nan"): 1},
            default=narwhals.col("col"),
        )
        .process(narwhals.from_native(df).lazy())
        .collect()
        .to_polars()
    )

    np.testing.assert_allclose(result["col"].to_numpy(), np.array([1, 1.5]))
