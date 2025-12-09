import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class StrLowercasePipe(DataProcessingPipe):
    target_column: str
    new_column_name: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.col(self.target_column)
            .str.to_lowercase()
            .alias(self.new_column_name)
        )
