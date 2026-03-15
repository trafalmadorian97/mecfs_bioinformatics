import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class TransposePipe(DataProcessingPipe):
    """
    Pipe to transpose a DataFrame
    Materializes data; Intended to run on small DataFrames.
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return narwhals.from_native(
            x.collect().to_pandas().transpose().reset_index()
        ).lazy()
