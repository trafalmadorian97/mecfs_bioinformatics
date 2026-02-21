from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class ShiftedLogPipe(DataProcessingPipe):
    base: int
    cols_to_exclude: Sequence[str]
    pseudocount: float = 1

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        for col_name in schema.keys():
            if col_name not in self.cols_to_exclude:
                x = x.with_columns(
                    (narwhals.col(col_name) + self.pseudocount)
                    .log(base=self.base)
                    .alias(col_name)
                )

        return x


@frozen
class ShiftedLogPipeInclude(DataProcessingPipe):
    base: int
    cols_to_include: Sequence[str]
    pseudocount: float = 1

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        for col_name in self.cols_to_include:
            x = x.with_columns(
                (narwhals.col(col_name) + self.pseudocount)
                .log(base=self.base)
                .alias(col_name)
            )
        return x
