import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class MoveColToFrontPipe(DataProcessingPipe):
    target_col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        cols = set(x.collect_schema().keys()) - set([self.target_col])
        cols = [self.target_col] + sorted(cols)
        return x.select(*cols)
