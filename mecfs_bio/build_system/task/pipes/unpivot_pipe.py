"""
Pipe to convert dataframe from wide to long format.
"""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class UnPivotPipe(DataProcessingPipe):
    on: str | list[str] | None
    index: str | list[str] | None
    variable_name: str
    value_name: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.unpivot(
            on=self.on,
            index=self.index,
            variable_name=self.variable_name,
            value_name=self.value_name,
        )
