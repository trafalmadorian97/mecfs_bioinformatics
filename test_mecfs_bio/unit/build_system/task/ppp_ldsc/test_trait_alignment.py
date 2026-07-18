"""Trait-to-context alignment: allele orientation, N sourcing, and the min-overlap guard."""

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.task.ppp_ldsc.trait_alignment import align_trait_to_context
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)


def _context(n: int = 10) -> pl.DataFrame:
    return pl.DataFrame(
        {
            GWASLAB_RSID_COL: [f"rs{i}" for i in range(n)],
            GWASLAB_EFFECT_ALLELE_COL: ["A"] * n,
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * n,
        }
    )


def test_orientation_match_flip_mismatch_and_absent():
    ctx = _context(4)
    trait = pl.DataFrame(
        {
            GWASLAB_RSID_COL: ["rs0", "rs1", "rs2"],  # rs3 absent from the trait
            GWASLAB_EFFECT_ALLELE_COL: ["A", "G", "C"],  # match, flipped, mismatch
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "A", "T"],
            GWASLAB_BETA_COL: [2.0, 2.0, 2.0],
            GWASLAB_SE_COL: [1.0, 1.0, 1.0],
            GWASLAB_SAMPLE_SIZE_COLUMN: [100.0, 100.0, 100.0],
        }
    )
    aligned = align_trait_to_context(
        trait, ctx, trait_total_sample_size=None, min_trait_snps=1
    )
    # rs0 same orientation -> +2; rs1 swapped -> -2; rs2 allele mismatch -> NaN; rs3 absent -> NaN.
    assert aligned.z[0] == pytest.approx(2.0)
    assert aligned.z[1] == pytest.approx(-2.0)
    assert np.isnan(aligned.z[2])
    assert np.isnan(aligned.z[3])
    assert aligned.n[0] == pytest.approx(100.0)
    assert np.isnan(aligned.n[2]) and np.isnan(aligned.n[3])


def test_constant_n_when_column_absent():
    ctx = _context(3)
    trait = pl.DataFrame(
        {
            GWASLAB_RSID_COL: ["rs0", "rs1", "rs2"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "A", "A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "G", "G"],
            GWASLAB_BETA_COL: [1.0, 1.0, 1.0],
            GWASLAB_SE_COL: [1.0, 1.0, 1.0],
        }
    )
    aligned = align_trait_to_context(
        trait, ctx, trait_total_sample_size=50_000, min_trait_snps=1
    )
    assert np.all(aligned.n == 50_000.0)


def test_missing_n_source_raises():
    ctx = _context(3)
    trait = pl.DataFrame(
        {
            GWASLAB_RSID_COL: ["rs0", "rs1", "rs2"],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "A", "A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "G", "G"],
            GWASLAB_BETA_COL: [1.0, 1.0, 1.0],
            GWASLAB_SE_COL: [1.0, 1.0, 1.0],
        }
    )
    with pytest.raises(AssertionError, match="trait_total_sample_size"):
        align_trait_to_context(
            trait, ctx, trait_total_sample_size=None, min_trait_snps=1
        )


def test_min_overlap_guard_raises():
    ctx = _context(10)
    trait = pl.DataFrame(
        {
            GWASLAB_RSID_COL: ["rs0"],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"],
            GWASLAB_BETA_COL: [1.0],
            GWASLAB_SE_COL: [1.0],
            GWASLAB_SAMPLE_SIZE_COLUMN: [100.0],
        }
    )
    with pytest.raises(AssertionError, match="harmonization failure"):
        align_trait_to_context(
            trait, ctx, trait_total_sample_size=None, min_trait_snps=5
        )
