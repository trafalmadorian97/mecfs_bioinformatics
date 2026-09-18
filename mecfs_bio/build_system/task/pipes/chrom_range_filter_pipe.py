"""Keep only rows in a single closed genomic interval. Mirrors the chrom-range filter that
HarmonizeGWASWithReferenceViaAlleles applied inline."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class ChromRangeFilterPipe(DataProcessingPipe):
    chrom: int
    start: int
    end: int
    chrom_col: str
    pos_col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(
            (narwhals.col(self.chrom_col) == self.chrom)
            & (narwhals.col(self.pos_col) >= self.start)
            & (narwhals.col(self.pos_col) <= self.end)
        )
