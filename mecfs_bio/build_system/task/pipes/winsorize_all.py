from typing import Sequence

import narwhals
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class WinsorizeAllPipe(DataProcessingPipe):
    max_value: float
    cols_to_exclude: Sequence[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        x_pl: pl.LazyFrame = x.collect().to_polars().lazy()
        schema = x_pl.collect_schema()
        for col_name in schema.keys():
            if col_name not in self.cols_to_exclude:
                x_pl = x_pl.with_columns(
                    pl.col(col_name).clip(upper_bound=self.max_value).alias(col_name)
                )

        return narwhals.from_native(x_pl)
