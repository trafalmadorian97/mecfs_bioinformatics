from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class CompositePipe(DataProcessingPipe):
    pipes: Sequence[DataProcessingPipe]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        for pipe in self.pipes:
            x = pipe.process(x)
        return x
