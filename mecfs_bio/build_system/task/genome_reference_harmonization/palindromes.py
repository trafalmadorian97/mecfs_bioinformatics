"""
Strand resolution of palindromic SNVs in untrusted tables, following gwaslab's rule.

Rows arrive oriented so that NEA is the reference base, so a panel record with
REF = NEA and ALT = EA describes the same variant. A variant is resolvable only when:

- its EAF is at most palindrome_maf_threshold, or at least 1 - palindrome_maf_threshold;
- the panel record exists;
- the panel MAF is at most panel_maf_threshold.

It is kept when EAF and panel AF lie on the same side of 0.5, and strand-flipped
(statistics flipped, alleles unchanged) otherwise.
"""

import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    PALINDROME_KEEP,
    PALINDROME_STRAND_FLIP,
    PALINDROME_UNRESOLVED,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

FREQUENCY_EPSILON = 1e-6
_PANEL_AF = "_panel_af"
_DECISION = "decision"
_KEYS = [GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]


def decide_palindrome_strands(
    rows: pl.DataFrame,
    panel: pl.DataFrame,
    options: GenomeReferenceHarmonizationOptions,
) -> pl.Series:
    """One PalindromeDecision per row of rows (POS, EA, NEA, EAF), in row order."""
    records = panel.select(
        GWASLAB_POS_COL,
        pl.col(PANEL_REF_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.col(PANEL_ALT_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(PANEL_AF_COL).cast(pl.Float64).alias(_PANEL_AF),
    )
    joined = rows.join(records, on=_KEYS, how="left", maintain_order="left")
    assert joined.height == rows.height, "duplicate panel records for a palindromic SNV"
    eaf = pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    af = pl.col(_PANEL_AF)
    threshold = options.palindrome_maf_threshold
    eaf_informative = (eaf <= threshold + FREQUENCY_EPSILON) | (
        eaf >= 1 - threshold - FREQUENCY_EPSILON
    )
    panel_informative = (
        pl.min_horizontal(af, 1 - af) <= options.panel_maf_threshold + FREQUENCY_EPSILON
    )
    resolvable = (eaf_informative & panel_informative).fill_null(False)
    same_side = ((af < 0.5) & (eaf < 0.5)) | ((af > 0.5) & (eaf > 0.5))
    return joined.select(
        pl.when(~resolvable)
        .then(pl.lit(PALINDROME_UNRESOLVED))
        .when(same_side)
        .then(pl.lit(PALINDROME_KEEP))
        .otherwise(pl.lit(PALINDROME_STRAND_FLIP))
        .alias(_DECISION)
    )[_DECISION]
