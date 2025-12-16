import narwhals
import pandas as pd
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.meta.read_spec.read_dataframe import ValidBackend
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class JointWithMemTablePipe(DataProcessingPipe):
    mem_table: pd.DataFrame
    keys_left: list[str]
    keys_right: list[str]
    backend: ValidBackend = "polars"

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        if self.backend == "polars":
            mem_table = narwhals.from_native(
                pl.from_pandas(
                    self.mem_table,
                )
            )
        else:
            raise ValueError("bad backend")

        return x.join(
            mem_table.lazy(), left_on=self.keys_left, right_on=self.keys_right
        )
