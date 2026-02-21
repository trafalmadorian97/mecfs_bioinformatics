import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class PivotPipe(DataProcessingPipe):
    """
    Transform a dataframe from long to wide format.
    """

    on: str | list[str]
    index: str | list[str] | None
    values: str | list[str] | None

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return (
            x.collect()
            .pivot(
                on=self.on,
                index=self.index,
                values=self.values,
            )
            .lazy()
        )
