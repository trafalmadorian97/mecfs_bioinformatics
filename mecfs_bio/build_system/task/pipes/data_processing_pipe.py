from abc import ABC, abstractmethod

import narwhals
import pandas as pd
import polars as pl


class DataProcessingPipe(ABC):
    """
    Abstract class representing a transformation of a lazy dataframe.
    """

    @abstractmethod
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        pass

    def process_pandas(self, x: pd.DataFrame) -> pd.DataFrame:
        lx = narwhals.from_native(x).lazy()
        return self.process(lx).collect().to_pandas()

    def process_eager_polars(self, x: pl.DataFrame) -> pl.DataFrame:
        lx = narwhals.from_native(x).lazy()
        return self.process(lx).collect().to_polars()
