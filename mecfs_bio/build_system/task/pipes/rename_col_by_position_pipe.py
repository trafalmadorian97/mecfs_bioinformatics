import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class RenameColByPositionPipe(DataProcessingPipe):
    col_position: int
    col_new_name: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        old_col = x.columns[self.col_position]
        return x.rename({old_col: self.col_new_name})
