import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class StrReplacePipe(DataProcessingPipe):
    target_column: str
    new_column_name: str
    replace_what: str
    replace_with: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.col(self.target_column)
            .str.replace_all(self.replace_what, self.replace_with)
            .alias(self.new_column_name)
        )
