"""Drop variants whose effect or non-effect allele is longer than one base."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
)


@frozen(slots=True)
class DropIndelsPipe(DataProcessingPipe):
    """
    Keep only variants whose EA and NEA are both a single base.

    Useful ahead of genome-reference harmonization for datasets with very long
    structural-variant alleles that downstream SNP-based analyses do not use.
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(
            _is_single_base(x, GWASLAB_EFFECT_ALLELE_COL)
            & _is_single_base(x, GWASLAB_NON_EFFECT_ALLELE_COL)
        )


def _is_single_base(x: narwhals.LazyFrame, col: str) -> narwhals.Expr:
    """
    Expression that is true where the allele in col has length one.

    A categorical column is tested against its single-base categories rather than cast to
    string: on a pandas-backed frame (as gwaslab produces) the cast allocates a fixed-width
    array of rows x (longest category) characters, which one long structural-variant allele
    makes too large to fit in memory.
    """
    if x.collect_schema()[col] != narwhals.Categorical:
        return narwhals.col(col).cast(narwhals.String).str.len_chars() == 1
    categories = x.select(col).collect().get_column(col).cat.get_categories()
    single_base = categories.filter(categories.str.len_chars() == 1).to_list()
    return narwhals.col(col).is_in(single_base)
