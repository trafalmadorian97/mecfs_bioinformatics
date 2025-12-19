import duckdb
import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class DuckdbMemLimitPipe(DataProcessingPipe):
    """
    Side-effecting pipe that sets the duckdb memory limit.
    Useful to apply prior to operations that consume a lot of memory.
    In these cases, counterintuitive, it often optimal to actually reduce the memory limit,
    as this will cause duckdb to more aggressively conserve memory.
    """

    limit_gb: int

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        duckdb.sql(f"SET memory_limit='{self.limit_gb}GB'")
        return x
