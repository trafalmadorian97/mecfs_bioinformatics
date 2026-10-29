"""
Colocalization of two traits at one locus by approximate Bayes factors (coloc.abf).

Follows Giambartolomei et al. 2014 and the coloc R package (approx.bf.estimates, combine.abf,
sdY.est). The method assumes at most one causal variant per trait in the locus and needs no LD.

Glossary:
- M: number of variants in the locus (both traits measured on the same M variants, same allele
  orientation).
- z = beta / se, V = se^2.
- W: prior standard deviation of the true effect. W = 0.15 * sdY for a quantitative trait (sdY is
  the trait's phenotypic standard deviation) and W = 0.2 for a case-control trait (log odds scale).
- r = W^2 / (W^2 + V).
- lABF: log approximate Bayes factor of a variant being causal versus no association,
  lABF = 0.5 * (log(1 - r) + r * z^2).

Hypotheses and their unnormalized log posteriors, with L1 = logsumexp(lABF1),
L2 = logsumexp(lABF2) and L12 = logsumexp(lABF1 + lABF2):
- H0, no association with either trait: 0.
- H1, trait 1 only: log p1 + L1.
- H2, trait 2 only: log p2 + L2.
- H3, both traits, distinct causal variants: log p1 + log p2 + log(exp(L1 + L2) - exp(L12)).
- H4, both traits, one shared causal variant: log p12 + L12.
p1, p2 and p12 are the prior probabilities that a given variant is causal for trait 1 only,
trait 2 only, or both.

When sdY of a quantitative trait is unknown it is estimated from the data: for a variant with
minor allele frequency f, sample size n and standard error se, se^2 is approximately
sdY^2 / (2 * n * f * (1 - f)). So sdY^2 is the no-intercept least squares coefficient of
2 * n * f * (1 - f) on 1 / se^2.
"""

from __future__ import annotations

from typing import Protocol

import numpy as np
from attrs import frozen
from scipy.special import logsumexp

_QUANT_PRIOR_SD_FACTOR = 0.15
_CASE_CONTROL_PRIOR_SD = 0.2
_SUM_TOLERANCE = 1e-9


@frozen(slots=True)
class QuantitativeTraitModel:
    """A quantitative trait; sd_y None means estimate sdY from allele frequency and N."""

    sd_y: float | None

    def __attrs_post_init__(self) -> None:
        assert self.sd_y is None or self.sd_y > 0, f"sd_y must be positive: {self.sd_y}"


@frozen(slots=True)
class CaseControlTraitModel:
    """A case-control trait with effects on the log odds scale."""

    case_fraction: float

    def __attrs_post_init__(self) -> None:
        assert 0 < self.case_fraction < 1, (
            f"case_fraction must be in (0, 1): {self.case_fraction}"
        )


TraitModel = QuantitativeTraitModel | CaseControlTraitModel


@frozen(slots=True)
class ColocPriors:
    p1: float = 1e-4
    p2: float = 1e-4
    p12: float = 1e-5


def _assert_vector(name: str, arr: np.ndarray, length: int) -> None:
    assert arr.ndim == 1, f"{name} must be 1-D, got shape {arr.shape}"
    assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"
    assert len(arr) == length, f"{name} has length {len(arr)}, expected {length}"


@frozen(slots=True)
class ColocInput:
    """
    Summary statistics of two traits on the same M variants, in the same allele orientation.

    maf and n of a trait are needed only when its model is quantitative with sd_y None.
    """

    beta1: np.ndarray
    se1: np.ndarray
    beta2: np.ndarray
    se2: np.ndarray
    maf1: np.ndarray | None
    n1: np.ndarray | None
    maf2: np.ndarray | None
    n2: np.ndarray | None
    model1: TraitModel
    model2: TraitModel

    def __attrs_post_init__(self) -> None:
        m = len(self.beta1)
        assert m > 0, "ColocInput needs at least one variant"
        for name, arr in [
            ("beta1", self.beta1),
            ("se1", self.se1),
            ("beta2", self.beta2),
            ("se2", self.se2),
        ]:
            _assert_vector(name, arr, m)
        assert np.all(self.se1 > 0) and np.all(self.se2 > 0), "se must be positive"
        for label, maf, n, model in [
            ("1", self.maf1, self.n1, self.model1),
            ("2", self.maf2, self.n2, self.model2),
        ]:
            if maf is not None:
                _assert_vector(f"maf{label}", maf, m)
            if n is not None:
                _assert_vector(f"n{label}", n, m)
            if isinstance(model, QuantitativeTraitModel) and model.sd_y is None:
                assert maf is not None and n is not None, (
                    f"trait {label} has no sd_y, so maf{label} and n{label} are required"
                )


@frozen(slots=True)
class ColocPosteriors:
    pp_h0: float
    pp_h1: float
    pp_h2: float
    pp_h3: float
    pp_h4: float
    n_variants: int

    def __attrs_post_init__(self) -> None:
        total = self.pp_h0 + self.pp_h1 + self.pp_h2 + self.pp_h3 + self.pp_h4
        assert abs(total - 1) < _SUM_TOLERANCE, f"posteriors sum to {total}"


class ColocMethod(Protocol):
    @property
    def name(self) -> str: ...

    def run(self, data: ColocInput) -> ColocPosteriors: ...


def estimate_sd_y(se: np.ndarray, maf: np.ndarray, n: np.ndarray) -> float:
    """Estimate sdY by no-intercept least squares of 2 * n * f * (1 - f) on 1 / se^2."""
    x = 1 / se**2
    y = 2 * n * maf * (1 - maf)
    keep = np.isfinite(x) & np.isfinite(y)
    assert keep.any(), "no finite variants to estimate sd_y from"
    coefficient = float(np.dot(x[keep], y[keep]) / np.dot(x[keep], x[keep]))
    return float(np.sqrt(coefficient))


def _prior_sd(
    model: TraitModel, se: np.ndarray, maf: np.ndarray | None, n: np.ndarray | None
) -> float:
    if isinstance(model, CaseControlTraitModel):
        return _CASE_CONTROL_PRIOR_SD
    sd_y = model.sd_y
    if sd_y is None:
        assert maf is not None and n is not None
        sd_y = estimate_sd_y(se=se, maf=maf, n=n)
    return _QUANT_PRIOR_SD_FACTOR * sd_y


def log_abf(beta: np.ndarray, se: np.ndarray, prior_sd: float) -> np.ndarray:
    """Per-variant lABF = 0.5 * (log(1 - r) + r * z^2), r = W^2 / (W^2 + V)."""
    v = se**2
    r = prior_sd**2 / (prior_sd**2 + v)
    z = beta / se
    return 0.5 * (np.log1p(-r) + r * z**2)


def _log_diff_exp(a: float, b: float) -> float:
    """log(exp(a) - exp(b)) for a >= b; -inf when they are equal."""
    if b >= a:
        return -np.inf
    return a + float(np.log1p(-np.exp(b - a)))


def combine_abf(
    labf1: np.ndarray, labf2: np.ndarray, priors: ColocPriors
) -> ColocPosteriors:
    l1 = float(logsumexp(labf1))
    l2 = float(logsumexp(labf2))
    l12 = float(logsumexp(labf1 + labf2))
    log_h = np.array(
        [
            0.0,
            np.log(priors.p1) + l1,
            np.log(priors.p2) + l2,
            np.log(priors.p1) + np.log(priors.p2) + _log_diff_exp(l1 + l2, l12),
            np.log(priors.p12) + l12,
        ]
    )
    pp = np.exp(log_h - logsumexp(log_h))
    return ColocPosteriors(
        pp_h0=float(pp[0]),
        pp_h1=float(pp[1]),
        pp_h2=float(pp[2]),
        pp_h3=float(pp[3]),
        pp_h4=float(pp[4]),
        n_variants=len(labf1),
    )


@frozen(slots=True)
class ColocAbf:
    """coloc.abf: single causal variant per trait, no LD."""

    priors: ColocPriors = ColocPriors()
    name: str = "coloc_abf"

    def run(self, data: ColocInput) -> ColocPosteriors:
        w1 = _prior_sd(data.model1, data.se1, data.maf1, data.n1)
        w2 = _prior_sd(data.model2, data.se2, data.maf2, data.n2)
        return combine_abf(
            labf1=log_abf(data.beta1, data.se1, w1),
            labf2=log_abf(data.beta2, data.se2, w2),
            priors=self.priors,
        )
