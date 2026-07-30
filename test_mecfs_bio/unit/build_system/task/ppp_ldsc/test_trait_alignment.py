"""Trait-to-context alignment: allele orientation, N sourcing, and the min-overlap guard."""

import narwhals
import numpy as np
import pandas as pd
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
        narwhals.from_native(trait.lazy()),
        ctx,
        trait_total_sample_size=None,
        min_trait_snps=1,
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
        narwhals.from_native(trait.lazy()),
        ctx,
        trait_total_sample_size=50_000,
        min_trait_snps=1,
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
    with pytest.raises(AssertionError):
        align_trait_to_context(
            narwhals.from_native(trait.lazy()),
            ctx,
            trait_total_sample_size=None,
            min_trait_snps=1,
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
    with pytest.raises(AssertionError):
        align_trait_to_context(
            narwhals.from_native(trait.lazy()),
            ctx,
            trait_total_sample_size=None,
            min_trait_snps=5,
        )


def _duplicated_trait(order: list[int]) -> pl.DataFrame:
    """A trait carrying rs0 twice with different effect sizes, emitted in the given row order."""
    rows = [
        {
            GWASLAB_RSID_COL: "rs0",
            GWASLAB_EFFECT_ALLELE_COL: "A",
            GWASLAB_NON_EFFECT_ALLELE_COL: "G",
            GWASLAB_BETA_COL: 1.0,
            GWASLAB_SE_COL: 1.0,
            GWASLAB_SAMPLE_SIZE_COLUMN: 100.0,
        },
        {
            GWASLAB_RSID_COL: "rs0",
            GWASLAB_EFFECT_ALLELE_COL: "A",
            GWASLAB_NON_EFFECT_ALLELE_COL: "G",
            GWASLAB_BETA_COL: 5.0,
            GWASLAB_SE_COL: 1.0,
            GWASLAB_SAMPLE_SIZE_COLUMN: 100.0,
        },
    ]
    return pl.DataFrame([rows[i] for i in order])


def test_duplicate_rsids_resolve_independently_of_input_order():
    """Which duplicate survives must depend on the data, not on the order the rows arrive in:
    a streaming scan makes no promise about that order, and an order-dependent answer would make
    the whole result -- and the asset built from it -- irreproducible."""
    ctx = _context(1)
    first = align_trait_to_context(
        narwhals.from_native(_duplicated_trait([0, 1]).lazy()),
        ctx,
        trait_total_sample_size=None,
        min_trait_snps=0,
    )
    reversed_input = align_trait_to_context(
        narwhals.from_native(_duplicated_trait([1, 0]).lazy()),
        ctx,
        trait_total_sample_size=None,
        min_trait_snps=0,
    )
    assert first.z[0] == reversed_input.z[0]
    # The smallest z-score wins, being first in sorted order.
    assert first.z[0] == pytest.approx(1.0)


def test_alignment_works_on_a_non_polars_backend():
    """The trait reaches alignment through a DataProcessingPipe, which may leave it on any
    narwhals backend, so the function must not assume the frame is polars-backed."""
    ctx = _context(2)
    trait = {
        GWASLAB_RSID_COL: ["rs0", "rs1"],
        GWASLAB_EFFECT_ALLELE_COL: ["A", "G"],  # match, flipped
        GWASLAB_NON_EFFECT_ALLELE_COL: ["G", "A"],
        GWASLAB_BETA_COL: [2.0, 2.0],
        GWASLAB_SE_COL: [1.0, 1.0],
        GWASLAB_SAMPLE_SIZE_COLUMN: [100.0, 100.0],
    }
    from_pandas = align_trait_to_context(
        narwhals.from_native(pd.DataFrame(trait)).lazy(),
        ctx,
        trait_total_sample_size=None,
        min_trait_snps=1,
    )
    from_polars = align_trait_to_context(
        narwhals.from_native(pl.DataFrame(trait).lazy()),
        ctx,
        trait_total_sample_size=None,
        min_trait_snps=1,
    )
    assert from_pandas.z == pytest.approx([2.0, -2.0])
    assert from_pandas.z == pytest.approx(from_polars.z)
    assert from_pandas.n == pytest.approx(from_polars.n)
