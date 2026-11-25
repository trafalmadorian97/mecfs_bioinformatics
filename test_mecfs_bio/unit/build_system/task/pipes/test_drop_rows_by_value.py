import polars as pl

from mecfs_bio.build_system.task.pipes.drop_rows_by_value import DropRowsByValue


def test_rows_with_listed_values_are_dropped() -> None:
    frame = pl.DataFrame({"kind": ["keep", "drop", "other", None]})
    result = DropRowsByValue(
        target_column="kind", values_to_drop=["drop"]
    ).process_eager_polars(frame)
    assert result["kind"].to_list() == ["keep", "other", None]
