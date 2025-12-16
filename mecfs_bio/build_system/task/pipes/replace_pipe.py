from typing import Mapping

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class ReplaceStrictPipe(DataProcessingPipe):
    target_column: str
    new_column_name: str
    replace_mapping: Mapping

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.col(self.target_column)
            .replace_strict(self.replace_mapping)
            .alias(self.new_column_name)
        )
