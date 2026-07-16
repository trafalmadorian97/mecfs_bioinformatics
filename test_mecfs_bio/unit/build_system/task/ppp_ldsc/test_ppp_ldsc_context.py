import numpy as np
import polars as pl

from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import (
    build_cis_mask,
    build_ppp_ldsc_context,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)
from mecfs_bio.constants.ppp_index_constants import PPP_INDEX_IS_STRAND_AMBIGUOUS_COL

# A row inside the extended MHC on hg38 (chr6, ~28 Mb) to prove MHC exclusion.
_MHC_POS = 28_000_000


def _index() -> pl.DataFrame:
    # Deliberately NOT position-sorted, so we can prove the context sorts by (CHR, POS)
    # and that row_pos refers back to original index rows.
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [2, 1, 1, 6, 1],
            GWASLAB_POS_COL: [300, 100, 200, _MHC_POS, 250],
            GWASLAB_RSID_COL: ["rs2", "rs1", "rs3", "rs_mhc", "rs_ambig"],
            PPP_INDEX_IS_STRAND_AMBIGUOUS_COL: [False, False, False, False, True],
        }
    )


def _ld() -> pl.DataFrame:
    # rs4 is absent from the index (inner join drops it); every index rsID is present.
    return pl.DataFrame(
        {
            "SNP": ["rs1", "rs2", "rs3", "rs_mhc", "rs_ambig", "rs4"],
            "L2": [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        }
    )


def test_context_sorts_filters_and_maps_row_positions():
    ctx = build_ppp_ldsc_context(
        _index(), _ld(), m_total=1000.0, drop_strand_ambiguous=True, exclude_mhc=True
    )
    # rs_ambig (strand-ambiguous) and rs_mhc (MHC) dropped -> rs1, rs3, rs2 remain,
    # sorted by (CHR, POS): chr1:100 (rs1), chr1:200 (rs3), chr2:300 (rs2).
    assert ctx.n_snps == 3
    assert list(ctx.chrom) == [1, 1, 2]
    assert list(ctx.pos) == [100, 200, 300]
    assert list(ctx.ld) == [1.0, 3.0, 2.0]
    assert list(ctx.wld) == [1.0, 3.0, 2.0]
    # row_pos points back to the original (unsorted) index rows: rs1=1, rs3=2, rs2=0.
    assert list(ctx.row_pos) == [1, 2, 0]
    assert ctx.m == 1000.0


def test_context_can_keep_strand_ambiguous_and_mhc():
    ctx = build_ppp_ldsc_context(
        _index(), _ld(), m_total=1.0, drop_strand_ambiguous=False, exclude_mhc=False
    )
    assert ctx.n_snps == 5  # nothing dropped (rs4 still absent from the index)


def test_build_cis_mask():
    ctx = build_ppp_ldsc_context(
        _index(), _ld(), m_total=1.0, drop_strand_ambiguous=True, exclude_mhc=True
    )
    # Gene on chr1 at [150, 160]; window 60 -> cis covers [90, 220]: rs1(100), rs3(200).
    mask = build_cis_mask(ctx, gene_chrom=1, gene_start=150, gene_end=160, window_bp=60)
    assert list(mask) == [True, True, False]
    # A gene on chr2 near the chr2 SNP (pos 300): window [220, 350] covers only rs2.
    other = build_cis_mask(
        ctx, gene_chrom=2, gene_start=280, gene_end=290, window_bp=60
    )
    assert list(other) == [False, False, True]
    # A chr1 SNP is never cis to the chr2 gene.
    assert list(
        build_cis_mask(ctx, gene_chrom=2, gene_start=150, gene_end=160, window_bp=60)
    ) == [False, False, False]
    assert isinstance(mask, np.ndarray) and mask.dtype == bool
