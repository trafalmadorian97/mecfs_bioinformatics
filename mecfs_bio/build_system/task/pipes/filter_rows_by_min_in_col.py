import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class FilterRowsByMinInCol(DataProcessingPipe):
    """
    Filer rows that do not exceed min value in  given col
    """

    min_value: float
    col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(
            narwhals.col(self.col) >= self.min_value,
        )
