from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class FilterRowsByValue(DataProcessingPipe):
    """
    Filer rows where column value lies in a given list
    """

    target_column: str
    valid_values: Sequence

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(
            narwhals.col(self.target_column).is_in(self.valid_values),
        )
