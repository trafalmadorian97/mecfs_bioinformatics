from typing import Sequence

import narwhals
from attrs import frozen
from narwhals.typing import UniqueKeepStrategy

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class UniquePipe(DataProcessingPipe):
    by: Sequence[str]
    keep: UniqueKeepStrategy
    order_by: str | Sequence[str] | None = None

    def __attrs_post_init__(self):
        if self.keep != "any":
            assert (
                self.order_by is not None
            )  # if we have a rule for deduplication based on ordering, we need to define an order

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.unique(subset=list(self.by), keep=self.keep, order_by=self.order_by)
