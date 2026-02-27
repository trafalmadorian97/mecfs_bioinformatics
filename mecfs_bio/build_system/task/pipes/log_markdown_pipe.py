import narwhals
import structlog

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe

logger = structlog.get_logger()


class LogMarkdownPipe(DataProcessingPipe):
    """
    For debugging
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        collected = x.collect().to_pandas()
        logger.debug(f"\n DataFrame:\n{collected.to_markdown(index=False)} \n")
        return x
