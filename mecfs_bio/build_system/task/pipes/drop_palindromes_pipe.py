"""Drop strand-ambiguous (palindromic) SNVs. Mirrors the palindrome handling that
HarmonizeGWASWithReferenceViaAlleles did inline, extracted so fine-mapping can apply it as a
gwas pipe once that task is no longer used. See is_palindromic_expr for the polars equivalent."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe

_PALINDROME_PAIRS = (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"))


@frozen
class DropPalindromesPipe(DataProcessingPipe):
    ea_col: str
    nea_col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        is_palindromic = narwhals.any_horizontal(
            (
                (narwhals.col(self.ea_col) == ea) & (narwhals.col(self.nea_col) == nea)
                for ea, nea in _PALINDROME_PAIRS
            ),
            ignore_nulls=True,
        )
        return x.filter(~is_palindromic)
