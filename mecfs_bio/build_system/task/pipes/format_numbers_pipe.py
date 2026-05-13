"""
Pipe to format a column of floating point values
"""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class FormatFloatNumbersPipe(DataProcessingPipe):
    """
    Pipe to set the format for floating point values in a given column.
    """

    format_str: str
    col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        collected = x.collect().to_pandas()
        result_col = []
        for value in collected[self.col].tolist():
            if _convertible_to_float(value):
                result_col.append(f"{float(value):{self.format_str}}")
            else:
                result_col.append(value)
        collected[self.col] = result_col
        return narwhals.from_native(collected).lazy()


def _convertible_to_float(val: str) -> bool:
    try:
        float(val)
        return True
    except (ValueError, TypeError):
        return False
