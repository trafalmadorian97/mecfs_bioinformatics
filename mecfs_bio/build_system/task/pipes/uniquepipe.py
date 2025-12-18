from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class UniquePipe(DataProcessingPipe):
    by: Sequence[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.unique(subset=list(self.by))
