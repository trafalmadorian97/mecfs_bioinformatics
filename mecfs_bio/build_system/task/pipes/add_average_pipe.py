from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class AddAveragePipe(DataProcessingPipe):
    cols_to_exclude: Sequence[str]
    avg_name: str = "Average"

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.mean_horizontal(narwhals.exclude(self.cols_to_exclude)).alias(
                self.avg_name
            )
        )
