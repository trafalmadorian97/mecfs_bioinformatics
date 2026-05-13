import narwhals

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


class ToPolarsPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return narwhals.from_native(x.collect().to_polars().lazy())
