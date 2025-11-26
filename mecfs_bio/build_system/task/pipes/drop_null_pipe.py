import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class DropNullsPipe(DataProcessingPipe):
    subset: None | list[str] = None

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        y = x.drop_nulls(subset=self.subset)
        return y
