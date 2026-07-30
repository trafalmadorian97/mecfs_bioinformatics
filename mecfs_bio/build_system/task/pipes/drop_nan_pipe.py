import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class DropNanPipe(DataProcessingPipe):
    cols: list[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        for col in self.cols:
            x = x.filter(~narwhals.col(col).is_nan())
        return x
