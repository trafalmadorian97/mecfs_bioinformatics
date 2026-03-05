import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_MLOG10P_COL, GWASLAB_P_COL


@frozen
class ComputeMlog10pIfNeededPipe(DataProcessingPipe):
    """
    Pipe to compute -log10(p).  There may be a narwhals bug when the underlying implementation is pandas
    """
    min_p_value: float = 1e-250

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        schema = x.collect_schema()
        if GWASLAB_MLOG10P_COL in schema:
            return x
        assert GWASLAB_P_COL in schema
        x = narwhals.from_native(x.collect().to_polars().lazy())  # workaround to fix issues with pandas backend
        result = x.with_columns((  -1* (narwhals.max_horizontal(  narwhals.col(GWASLAB_P_COL).cast(narwhals.dtypes.Float64()),narwhals.lit(self.min_p_value),  ).log(base=10,))).alias(GWASLAB_MLOG10P_COL))
        return result
