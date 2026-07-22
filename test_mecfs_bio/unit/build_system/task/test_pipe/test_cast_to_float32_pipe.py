import narwhals
import polars as pl

from mecfs_bio.build_system.task.pipes.cast_to_float32_pipe import (
    CastFloatsToFloat32Pipe,
)


def test_cast_to_float32():
    df = pl.DataFrame({"col": [1.4, 2.5]}, schema={"col": pl.Float64()})
    result = CastFloatsToFloat32Pipe().process(narwhals.from_native(df).lazy())
    assert result.collect_schema()["col"] == narwhals.dtypes.Float32()
