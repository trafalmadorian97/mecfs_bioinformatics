import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class CastIntsToInt32Pipe(DataProcessingPipe):
    """
    Convert all int64 to int32
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.selectors.by_dtype(narwhals.Int64).cast(narwhals.Int32)
        )
