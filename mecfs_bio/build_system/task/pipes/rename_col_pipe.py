import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class RenameColPipe(DataProcessingPipe):

    old_name:str
    new_name:str
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.rename({self.old_name:self.new_name})
