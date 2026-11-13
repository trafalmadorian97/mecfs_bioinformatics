from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class DropRowsByValue(DataProcessingPipe):
    """
    Drop rows whose column value lies in a given list. Rows with a null value are kept.
    """

    target_column: str
    values_to_drop: Sequence

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        col = narwhals.col(self.target_column)
        return x.filter(~col.is_in(self.values_to_drop) | col.is_null())
