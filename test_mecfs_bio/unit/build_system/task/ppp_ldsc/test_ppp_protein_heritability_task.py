"""Unit tests for the sample-size helper used when PppProteinHeritabilityTask reads N from
the slim parquet (sample_size_task=None) instead of a per-protein sample-size table."""

import math

import numpy as np
import pytest

from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_heritability_task import (
    constant_sample_size,
)


def test_constant_sample_size_ignores_nan_absent_variants():
    n_at_context = np.array([np.nan, 42000.0, 42000.0, np.nan, 42000.0])
    assert constant_sample_size(n_at_context, "PROTEIN") == pytest.approx(42000.0)


def test_constant_sample_size_raises_when_not_constant():
    n_at_context = np.array([42000.0, 42001.0, np.nan])
    with pytest.raises(AssertionError, match="distinct sample sizes"):
        constant_sample_size(n_at_context, "PROTEIN")


def test_constant_sample_size_raises_when_all_absent():
    # No finite value -> zero distinct sizes -> fail (cannot recover N).
    n_at_context = np.array([np.nan, np.nan])
    with pytest.raises(AssertionError, match="distinct sample sizes"):
        constant_sample_size(n_at_context, "PROTEIN")
    assert math.isnan(n_at_context[0])
