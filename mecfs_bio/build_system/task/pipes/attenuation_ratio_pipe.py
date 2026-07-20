import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class AttenuationRatioPipe(DataProcessingPipe):
    """
    Pipe to compute attenuation ratio.
    Useful for evaluating risk of population stratification
    """

    mean_chi_col: str
    intercept_col: str
    ratio_col: str = "ratio"

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.with_columns(
            narwhals.when(
                (narwhals.col(self.mean_chi_col) >= 1.02)
                & (narwhals.col(self.intercept_col) >= 1),
            )
            .then(
                (narwhals.col(self.intercept_col) - 1)
                / (narwhals.col(self.mean_chi_col) - 1)
            )
            .otherwise(float("nan"))
            .alias(self.ratio_col)
        )
