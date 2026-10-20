import narwhals
import polars as pl

from mecfs_bio.build_system.task.pipes.cast_to_float32_pipe import (
    CastFloatsToFloat32Pipe,
)


def test_cast_to_float32() -> None:
    df = pl.DataFrame({"col": [1.4, 2.5]}, schema={"col": pl.Float64()})
    result = CastFloatsToFloat32Pipe().process(narwhals.from_native(df).lazy())
    assert result.collect_schema()["col"] == narwhals.dtypes.Float32()


def test_cast_to_float32_respects_exclusions() -> None:
    df = pl.DataFrame(
        {"cast_me": [1.4, 2.5], "keep_me": [1.4, 2.5]},
        schema={"cast_me": pl.Float64(), "keep_me": pl.Float64()},
    )
    result = CastFloatsToFloat32Pipe(cols_to_exclude=["keep_me"]).process(
        narwhals.from_native(df).lazy()
    )
    schema = result.collect_schema()
    assert schema["cast_me"] == narwhals.dtypes.Float32()
    assert schema["keep_me"] == narwhals.dtypes.Float64()
