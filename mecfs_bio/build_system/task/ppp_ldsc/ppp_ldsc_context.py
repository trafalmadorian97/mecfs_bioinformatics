"""
Shared, protein-invariant context for LDSC over the UKB-PPP database.

Every per-protein slim file stores beta/se in the SAME variant-index row order, so the
alignment of index variants to LD scores, the regression LD scores (ld), the total
reference-SNP count M, and the genome-sorted order used for contiguous jackknife blocks are
all identical across proteins. We therefore build them ONCE and reuse them for every
protein (and, later, for cross-trait rg and LCV).

The context also carries each retained variant's chromosome and hg38 position so a per-
protein cis mask (variants within +/- a window of the protein's gene) can be computed, and
its row position in the slim files so a protein's beta/se can be gathered by positional
index.

Built by intersecting the variant index (rsID) with the reference LD scores (SNP), then
optionally dropping strand-ambiguous variants and the extended MHC region (on the index's
primary hg38 POS), and sorting by (chromosome, position).

The reference LD scores arrive as a consolidated dataframe (ConsolidateLDScoresTask), which also
carries the common-variant count M_5_50 the regression normalizes by, so this module reads
dataframes rather than the LDSC authors' per-chromosome file layout.
"""

from __future__ import annotations

import narwhals
import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_LD_SCORE_COL,
    LD_SCORE_RSID_COL,
    total_m_5_50,
)
from mecfs_bio.constants.genomic_coordinate_constants import extended_mhc_interval
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_database_constants import (
    PPP_INDEX_IS_STRAND_AMBIGUOUS_COL,
)

# The reference LD scores are hg19-coordinate; the index (and thus our exclusions) use hg38.
_ROW_POS_COL = "__ctx_row__"


@frozen
class PppLdscContext:
    """The shared regression SNP set, genome-sorted. All arrays are parallel (one entry per
    retained variant, in (chromosome, position) order).

    We do not carry a separate regression-weight LD score (wLD): with sep_weights=False the
    weight LD score equals ld, and the batched kernel reuses ld directly. Reinstate a wld
    array here only if a future path (e.g. externally supplied regression weights) needs it.
    """

    row_pos: (
        np.ndarray
    )  # (S,) int64: position of each variant in the slim-file row order
    ld: np.ndarray  # (S,) float: LD score (L2)
    chrom: np.ndarray  # (S,) int64: chromosome
    pos: np.ndarray  # (S,) int64: hg38 position (for cis masks)
    m: float  # total reference-SNP count (sum of M_5_50)

    def __attrs_post_init__(self) -> None:
        # Make invalid states unrepresentable: every array is 1-D, the same length, and the
        # expected dtype. Constructing a context that violates these is a bug, not a warning.
        s = self.row_pos.shape[0]
        for name, arr, kind in (
            ("row_pos", self.row_pos, "i"),
            ("ld", self.ld, "f"),
            ("chrom", self.chrom, "i"),
            ("pos", self.pos, "i"),
        ):
            assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
            assert arr.shape[0] == s, f"{name} has length {arr.shape[0]}, expected {s}"
            assert arr.dtype.kind == kind, (
                f"{name} must be {kind!r}-kind, got dtype {arr.dtype}"
            )

    @property
    def n_snps(self) -> int:
        return int(self.row_pos.shape[0])


def build_ppp_ldsc_context(
    index_df: pl.DataFrame,
    ld_df: pl.DataFrame,
    *,
    drop_strand_ambiguous: bool = True,
    exclude_mhc: bool = True,
) -> PppLdscContext:
    """Intersect the variant index with the LD scores and assemble the shared context.

    index_df: the variant index, with columns CHR, POS (hg38), rsID, and
        is_strand_ambiguous, in its canonical row order (row i aligns to row i of every
        slim file).
    ld_df: consolidated reference LD scores with columns CHR, SNP (rsID), L2 and M_5_50.

    M is taken from the LD scores BEFORE the join, since it counts reference variants and must
    not shrink to the subset the index happens to cover.
    """
    m_total = total_m_5_50(narwhals.from_native(ld_df, eager_only=True))
    joined = (
        index_df.with_row_index(_ROW_POS_COL)
        .join(
            ld_df.select(LD_SCORE_RSID_COL, LD_SCORE_LD_SCORE_COL),
            left_on=GWASLAB_RSID_COL,
            right_on=LD_SCORE_RSID_COL,
            how="inner",
        )
        .sort([GWASLAB_CHROM_COL, GWASLAB_POS_COL])
    )
    if drop_strand_ambiguous:
        joined = joined.filter(~pl.col(PPP_INDEX_IS_STRAND_AMBIGUOUS_COL))
    if exclude_mhc:
        mhc = extended_mhc_interval("38")
        in_mhc = (
            (pl.col(GWASLAB_CHROM_COL) == mhc.chrom)
            & (pl.col(GWASLAB_POS_COL) >= mhc.start)
            & (pl.col(GWASLAB_POS_COL) <= mhc.end)
        )
        joined = joined.filter(~in_mhc)

    return PppLdscContext(
        row_pos=joined[_ROW_POS_COL].to_numpy().astype(np.int64),
        ld=joined[LD_SCORE_LD_SCORE_COL].to_numpy().astype(float),
        chrom=joined[GWASLAB_CHROM_COL].to_numpy().astype(np.int64),
        pos=joined[GWASLAB_POS_COL].to_numpy().astype(np.int64),
        m=m_total,
    )


def build_cis_mask(
    context: PppLdscContext,
    gene_chrom: int,
    gene_start: int,
    gene_end: int,
    window_bp: int,
) -> np.ndarray:
    """Boolean over the context SNP set, True for CIS variants (to be excluded): those on
    the gene's chromosome within window_bp of the gene body [start, end]."""
    return (
        (context.chrom == gene_chrom)
        & (context.pos >= gene_start - window_bp)
        & (context.pos <= gene_end + window_bp)
    )
