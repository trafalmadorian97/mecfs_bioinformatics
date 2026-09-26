from typing import Sequence

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen(slots=True)
class FilterRowsByMin(DataProcessingPipe):
    """
    Filer rows that do not exceed min value in at least one column
    """

    min_value: float
    exclude_columns: Sequence[str]

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        x_pd = x.collect().to_pandas()
        cols = x_pd.columns
        cols = sorted(set(cols) - set(self.exclude_columns))
        x_pd = x_pd.loc[(x_pd[cols] >= self.min_value).any(axis=1), :]
        return narwhals.from_native(x_pd).lazy()
