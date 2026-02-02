import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class HeadPipe(DataProcessingPipe):
    """
    Pipe to extract first rows of a dataframe

    Can be useful in writing tests that operate on a subset of a large dataset.
    """

    num_rows: int

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.head(self.num_rows)
