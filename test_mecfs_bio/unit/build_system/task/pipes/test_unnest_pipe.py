import polars as pl

from mecfs_bio.build_system.task.pipes.unnest_pipe import UNNestPipe


def test_list_column_becomes_one_row_per_element() -> None:
    frame = pl.DataFrame({"gene": ["A", "B"], "ids": [["P1", "P2"], ["P3"]]})
    result = UNNestPipe(col_to_unnest="ids").process_eager_polars(frame)
    assert result.rows() == [("A", "P1"), ("A", "P2"), ("B", "P3")]
