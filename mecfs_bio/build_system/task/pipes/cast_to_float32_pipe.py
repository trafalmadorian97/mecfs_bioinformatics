from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class CastFloatsToFloat32Pipe(DataProcessingPipe):
    """
    Convert all float64 columns to float32, except those named in cols_to_exclude.
    """

    cols_to_exclude: Sequence[str] = ()

    def __attrs_post_init__(self) -> None:
        assert not isinstance(self.cols_to_exclude, str), (
            "Pass a list of column names to exclude, not a string"
        )

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        missing = set(self.cols_to_exclude) - set(schema.names())
        assert not missing, (
            f"Excluded columns not present in frame: {sorted(missing)}.  List of present columns: {sorted(schema.names())}"
        )
        cols_to_cast = [
            col_name
            for col_name, dtype in schema.items()
            if dtype == narwhals.Float64 and col_name not in self.cols_to_exclude
        ]
        return x.with_columns(narwhals.col(cols_to_cast).cast(narwhals.Float32))
