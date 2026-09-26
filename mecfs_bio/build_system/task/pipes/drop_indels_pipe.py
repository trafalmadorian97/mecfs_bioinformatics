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
            (
                narwhals.col(GWASLAB_EFFECT_ALLELE_COL)
                .cast(narwhals.String)
                .str.len_chars()
                == 1
            )
            & (
                narwhals.col(GWASLAB_NON_EFFECT_ALLELE_COL)
                .cast(narwhals.String)
                .str.len_chars()
                == 1
            )
        )
