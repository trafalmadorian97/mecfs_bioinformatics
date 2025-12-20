from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class SortPipe(DataProcessingPipe):
    by: Sequence[str]
    desc: Sequence[bool] | bool = False

    def __attrs_post_init__(self):
        if self.desc is not None:
            assert len(self.by) == len(self.desc)

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.sort(by=self.by, descending=self.desc)
