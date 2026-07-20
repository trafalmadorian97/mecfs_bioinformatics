import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class CastToFloat32Pipe(DataProcessingPipe):
    """
    Convert all float64 columns to float32
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.selectors.by_dtype(narwhals.Float64).cast(narwhals.Float32)
        )
