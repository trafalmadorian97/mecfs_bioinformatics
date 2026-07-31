import narwhals
import polars as pl
from attrs import frozen
from statsmodels.stats.multitest import multipletests

from mecfs_bio.build_system.task.multiple_testing_table_task import (
    REJECT_NULL_LABEL,
    Method,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class MultipleTestingPipe(DataProcessingPipe):
    p_col: str
    alpha: float = 0.05
    method: Method = "bonferroni"
    reject_name: str = REJECT_NULL_LABEL

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        df = x.collect().to_polars()
        p = df[self.p_col].to_numpy()
        reject, _, _, _ = multipletests(pvals=p, alpha=self.alpha, method=self.method)
        reject_ser = pl.Series(name=self.reject_name, values=reject)
        df = df.with_columns(reject_ser)
        return narwhals.from_native(df).lazy()
