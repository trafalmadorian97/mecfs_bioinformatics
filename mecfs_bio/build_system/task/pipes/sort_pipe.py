import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class SortPipe(DataProcessingPipe):
    by: list[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.sort(by=self.by)
