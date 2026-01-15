from abc import ABC, abstractmethod

import narwhals
import pandas as pd


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
