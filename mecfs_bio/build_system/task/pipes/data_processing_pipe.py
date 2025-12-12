from abc import ABC, abstractmethod

import narwhals


class DataProcessingPipe(ABC):
    """
    Abstract class representing a transformation of a lazy dataframe.
    """

    @abstractmethod
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        pass
