import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class SplitColPipe(DataProcessingPipe):
    col_to_split: str
    split_by: str
    new_col_name: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.col(self.col_to_split)
            .str.split(self.split_by)
            .alias(self.new_col_name)
        )
