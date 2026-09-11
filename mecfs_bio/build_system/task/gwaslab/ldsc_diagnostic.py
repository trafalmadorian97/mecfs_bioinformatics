"""
The arithmetic behind the LD-score-regression diagnostic plot.

LD score regression models the per-variant chi-square as a straight line in the LD score:
E[chi^2_j] = intercept + (N * h2 / M) * ld_j, where N is the sample size, M the number of
reference variants, and h2 the observed-scale heritability. The diagnostic plot bins variants
by LD score and shows each bin's mean chi-square against that line, so a bin that sits well off
the line, or an intercept that drifts away from one, is visible by eye.

This module holds only the pure math: binning variants into equal-count LD-score bins, the
fitted line, and the rescaling that turns the chi-square axis into heritability-slope units. It
takes plain numpy arrays and returns plain values, so it can be tested without reading any asset
or drawing anything.
"""

from __future__ import annotations

import numpy as np
from attrs import frozen


@frozen
class LdscDiagnosticBins:
    """Per-bin summaries of variants grouped into equal-count LD-score bins, ordered by
    ascending LD score. Every array is 1-D of length n_bins."""

    mean_ld: np.ndarray
    mean_chi2: np.ndarray
    se_chi2: np.ndarray
    count: np.ndarray

    def __attrs_post_init__(self) -> None:
        n = self.mean_ld.shape[0]
        for name, arr, kind in (
            ("mean_ld", self.mean_ld, "f"),
            ("mean_chi2", self.mean_chi2, "f"),
            ("se_chi2", self.se_chi2, "f"),
            ("count", self.count, "i"),
        ):
            assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
            assert arr.shape[0] == n, f"{name} has length {arr.shape[0]}, expected {n}"
            assert arr.dtype.kind == kind, (
                f"{name} must be {kind!r}-kind, got dtype {arr.dtype}"
            )


def bin_by_ld_score(
    chi2: np.ndarray, ld_score: np.ndarray, n_bins: int
) -> LdscDiagnosticBins:
    """Group variants into n_bins equal-count bins by ascending LD score and summarize each.

    Equal-count (rather than equal-width) bins put a comparable number of variants in every
    point, so each mean carries similar noise -- matching the binned plots in the original LDSC
    papers. Bin sizes differ by at most one when n_bins does not divide the variant count."""
    assert chi2.shape == ld_score.shape, (
        f"chi2 and ld_score must have the same shape, got {chi2.shape} and {ld_score.shape}"
    )
    assert n_bins >= 1, f"n_bins must be at least 1, got {n_bins}"

    order = np.argsort(ld_score, kind="stable")
    ld_sorted = ld_score[order]
    chi2_sorted = chi2[order]

    mean_ld = []
    mean_chi2 = []
    se_chi2 = []
    count = []
    for ld_bin, chi2_bin in zip(
        np.array_split(ld_sorted, n_bins), np.array_split(chi2_sorted, n_bins)
    ):
        c = chi2_bin.shape[0]
        mean_ld.append(float(ld_bin.mean()))
        mean_chi2.append(float(chi2_bin.mean()))
        # Sample std (ddof=1) needs at least two values; a singleton bin has no spread.
        se_chi2.append(float(chi2_bin.std(ddof=1) / np.sqrt(c)) if c > 1 else 0.0)
        count.append(c)

    return LdscDiagnosticBins(
        mean_ld=np.array(mean_ld, dtype=float),
        mean_chi2=np.array(mean_chi2, dtype=float),
        se_chi2=np.array(se_chi2, dtype=float),
        count=np.array(count, dtype=np.int64),
    )


def chi2_to_heritability_units(chi2: np.ndarray, n: float, m: float) -> np.ndarray:
    """Rescale chi-square to chi^2 * M / N, the units in which the fitted line's slope equals
    the observed-scale heritability. This is the secondary y-axis, and because M and N are
    divided out it is comparable across traits of different sample size."""
    return chi2 * m / n


def fit_line_chi2(
    intercept: float, h2_obs: float, n: float, m: float, ld: np.ndarray
) -> np.ndarray:
    """The LDSC fitted line in chi-square units: intercept + (N * h2 / M) * ld."""
    slope = n * h2_obs / m
    return intercept + slope * ld
