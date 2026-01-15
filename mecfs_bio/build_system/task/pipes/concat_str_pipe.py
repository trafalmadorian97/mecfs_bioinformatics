from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class ConcatStrPipe(DataProcessingPipe):
    target_cols: Sequence[str]
    sep: str
    new_col_name: str
    ignore_nulls: bool = True

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.concat_str(
                [narwhals.col(item) for item in self.target_cols],
                separator=self.sep,
                ignore_nulls=self.ignore_nulls,
            ).alias(self.new_col_name)
        )
