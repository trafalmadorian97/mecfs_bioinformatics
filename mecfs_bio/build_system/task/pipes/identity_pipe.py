import narwhals

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


class IdentityPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x
