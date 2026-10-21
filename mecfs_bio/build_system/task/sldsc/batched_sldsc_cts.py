"""
Batched cell-type-specific stratified LD score regression.

Reproduces ldsc's cell_type_specific for many models at once. ldsc fits each model with its Hsq
class (Hsq = "h-squared", ldsc's heritability estimator in gwaslab/extension/ldsc/ldsc_regressions.py,
which fits LD score regression with a block-jackknife standard error); cell_type_specific calls it
with old_weights=True, meaning a single weighted least squares with no iterative reweighting and
no two-step intercept estimate. Each model regresses chi-square on [target, controls, baseline]
LD scores scaled by N / mean(N), plus an intercept, and reports the target coefficient with a
200-block jackknife standard error.

M_5_50: for each annotation, the sum of its values over reference-panel SNPs with minor allele
frequency above 5% (for a binary annotation, the number of common SNPs it contains). In the LDSC
model E[chi2_j] = N_j * sum_c tau_c * l(j, c) + N_j * a + 1, where l(j, c) is the LD score of SNP j
with annotation c and tau_c is the per-SNP heritability contributed by annotation c, M_5_50
converts per-SNP effects into the heritability attributed to an annotation
(h2_c = tau_c * M_5_50_c). Here the total M_5_50 over a model's annotations scales the aggregate
heritability estimate that shapes the regression weights.

Models sharing a control set share the design apart from the target column, so they are
processed together as a control group. In exact mode each model keeps its own regression weights
(the weights depend on the model's total LD and total M_5_50, which include the target); the
per-block normal equations of a chunk of models are formed with a few large matrix products
rather than one regression at a time. ldsc normalizes the square-root weights to sum to one before
weighting the design and response; that scalar cancels in every solve, so the effective weight
here is the raw ldsc weight.

In shared_controls mode every model in a control group uses the weight of the model without a
target column, so the shared Gram matrix is formed once and only the target's border changes; the
target coefficient follows from the Frisch-Waugh-Lovell identity.

Dimension names used throughout this module:
    S: number of regression SNPs (rows of every per-SNP array).
    B: number of baseline LD-score columns.
    U: number of LD-score columns in the annotation set (targets and controls).
    K: number of models (one reported target column each).
    k: number of models in the chunk currently being processed (k <= chunk_size).
    q: number of shared design columns of a control group: its controls, the baseline columns,
       and the intercept (q = number of controls + B + 1).
    p: number of design columns of one model: its target plus the shared columns (p = q + 1).
    n_blocks: number of jackknife blocks, min(S, 200).
"""

from collections import defaultdict
from typing import Sequence

import numpy as np
from attrs import frozen
from scipy.linalg import lu_factor, lu_solve
from scipy.stats import norm

from mecfs_bio.constants.sldsc_constants import SLDSCWeighting

LDSC_DEFAULT_N_BLOCKS = 200
DEFAULT_CHUNK_SIZE = 64
_NULL_INTERCEPT = 1.0


@frozen(slots=True)
class SLDSCRegressionInputs:
    """
    Regression SNPs in ldsc row order (the jackknife blocks depend on it). Dimensions S, B, U
    are defined in the module docstring.

    chi2, n, w_ld: (S,) chi-square, integer-valued sample size, weight LD score.
    baseline_ld: (S, B) baseline LD scores; baseline_m: (B,) baseline M_5_50.
    annotation_ld: (S, U) target and control LD scores (float32 allowed); annotation_m: (U,)
        their M_5_50.
    """

    chi2: np.ndarray
    n: np.ndarray
    w_ld: np.ndarray
    baseline_ld: np.ndarray
    baseline_m: np.ndarray
    annotation_ld: np.ndarray
    annotation_m: np.ndarray

    def __attrs_post_init__(self) -> None:
        s = self.chi2.shape[0]
        for name, arr in (("chi2", self.chi2), ("n", self.n), ("w_ld", self.w_ld)):
            assert arr.ndim == 1 and arr.shape[0] == s, (
                f"{name} must have shape ({s},), got {arr.shape}"
            )
            assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"
        for name, ld, m in (
            ("baseline", self.baseline_ld, self.baseline_m),
            ("annotation", self.annotation_ld, self.annotation_m),
        ):
            assert ld.ndim == 2 and ld.shape[0] == s, (
                f"{name}_ld must have shape ({s}, columns), got {ld.shape}"
            )
            assert ld.dtype.kind == "f", f"{name}_ld must be float, got {ld.dtype}"
            assert m.ndim == 1 and m.shape[0] == ld.shape[1], (
                f"{name}_m must have shape ({ld.shape[1]},), got {m.shape}"
            )
        assert np.all(np.isfinite(self.chi2)), (
            "chi2 must be finite (filter before the kernel)"
        )


@frozen(slots=True)
class SLDSCModel:
    """target_index and control_indices are columns of SLDSCRegressionInputs.annotation_ld."""

    label: str
    target_index: int
    control_indices: tuple[int, ...]

    def __attrs_post_init__(self) -> None:
        assert self.target_index not in self.control_indices, (
            f"{self.label}: target among controls"
        )


@frozen(slots=True)
class BatchedSLDSCResult:
    """
    Per-model target coefficient, jackknife SE, and one-sided p-value, each (K,), in input model
    order (K = number of models).
    """

    labels: tuple[str, ...]
    coef: np.ndarray
    coef_se: np.ndarray
    p: np.ndarray

    def __attrs_post_init__(self) -> None:
        k = len(self.labels)
        for name, arr in (
            ("coef", self.coef),
            ("coef_se", self.coef_se),
            ("p", self.p),
        ):
            assert arr.ndim == 1 and arr.shape[0] == k, (
                f"{name} must have shape ({k},), got {arr.shape}"
            )
            assert arr.dtype.kind == "f", f"{name} must be float, got {arr.dtype}"


@frozen(slots=True)
class _SharedDesign:
    """
    Columns common to every model of one control group.

    ld_sum: (S,) row-sum of the unscaled control and baseline LD scores (S = regression SNPs).
    m_sum: total M_5_50 of the controls and baseline columns.
    design: (S, q) scaled design: controls and baseline LD scores times N / mean(N), then a final
        intercept column of ones (q = number of controls + number of baseline columns + 1).
    """

    ld_sum: np.ndarray
    m_sum: float
    design: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.ld_sum.ndim == 1, f"ld_sum must be 1-D, got {self.ld_sum.shape}"
        assert self.design.ndim == 2 and self.design.shape[0] == self.ld_sum.shape[0], (
            f"design must have shape ({self.ld_sum.shape[0]}, q), got {self.design.shape}"
        )
        assert self.ld_sum.dtype.kind == "f" and self.design.dtype.kind == "f", (
            "shared design arrays must be float"
        )


@frozen(slots=True)
class _BlockNormalEquations:
    """
    Per-block weighted normal equations of k models with p design columns each, target first.

    gram: (n_blocks, k, p, p) per-block X'WX.
    rhs: (n_blocks, k, p) per-block X'Wy.
    """

    gram: np.ndarray
    rhs: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.gram.ndim == 4 and self.gram.shape[2] == self.gram.shape[3], (
            f"gram must have shape (n_blocks, k, p, p), got {self.gram.shape}"
        )
        assert self.rhs.shape == self.gram.shape[:3], (
            f"rhs must have shape {self.gram.shape[:3]}, got {self.rhs.shape}"
        )


@frozen(slots=True)
class _TargetEstimates:
    """
    Target coefficient of k models, before division by mean(N).

    est: (k,) estimate on all regression SNPs.
    delete: (n_blocks, k) estimate with each jackknife block left out.
    """

    est: np.ndarray
    delete: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.est.ndim == 1, f"est must be 1-D, got {self.est.shape}"
        assert self.delete.ndim == 2 and self.delete.shape[1] == self.est.shape[0], (
            f"delete must have shape (n_blocks, {self.est.shape[0]}), got {self.delete.shape}"
        )


def ldsc_block_separators(n_snps: int, n_blocks: int) -> np.ndarray:
    """ldsc's Jackknife.get_separators: floor(linspace(0, n_snps, n_blocks + 1))."""
    return np.floor(np.linspace(0, n_snps, n_blocks + 1)).astype(int)


def ldsc_hsq_weights(
    ld_tot: np.ndarray,
    w_ld: np.ndarray,
    n: np.ndarray,
    m_tot: np.ndarray,
    chi2: np.ndarray,
) -> np.ndarray:
    """
    Regression weights of ldsc's Hsq estimator (old_weights=True path), for k models at once.

    ld_tot: (S, k) unscaled total LD score of each model (sum over its annotations).
    w_ld: (S,) weight LD score; n: (S,) sample size; chi2: (S,) chi-square.
    m_tot: (k,) total M_5_50 of each model.

    Returns (S, k) weights het * oc, where het = 1 / (2 (1 + h2_agg N max(ld_tot, 1) / m_tot)^2)
    corrects for heteroskedasticity, oc = 1 / max(w_ld, 1) for over-counting of correlated SNPs,
    and h2_agg = clip(m_tot (mean chi2 - 1) / mean(ld_tot N), 0, 1). ldsc weights the design and
    response by the square root of this, so it is the effective least-squares weight.
    """
    h2_agg = np.clip(
        m_tot * (chi2.mean() - _NULL_INTERCEPT) / (ld_tot * n[:, None]).mean(axis=0),
        0.0,
        1.0,
    )
    c = h2_agg[None, :] * n[:, None] / m_tot[None, :]
    het = 1.0 / (2.0 * (_NULL_INTERCEPT + c * np.maximum(ld_tot, 1.0)) ** 2)
    return het / np.maximum(w_ld, 1.0)[:, None]


def batched_sldsc_cts(
    inputs: SLDSCRegressionInputs,
    models: Sequence[SLDSCModel],
    *,
    weighting: SLDSCWeighting,
    n_blocks: int = LDSC_DEFAULT_N_BLOCKS,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
) -> BatchedSLDSCResult:
    assert models, "need at least one model"
    s = inputs.chi2.shape[0]
    n_blocks = min(s, n_blocks)
    separators = ldsc_block_separators(s, n_blocks)
    n_bar = float(inputs.n.mean())
    scale = inputs.n / n_bar
    est = np.empty(len(models))
    delete = np.empty((n_blocks, len(models)))
    fit_group = _exact_group if weighting == "exact" else _shared_controls_group
    for controls, positions in _group_by_controls(models).items():
        group = fit_group(
            inputs=inputs,
            shared=_shared_design(inputs=inputs, controls=controls, scale=scale),
            targets=[models[i].target_index for i in positions],
            scale=scale,
            separators=separators,
            chunk_size=chunk_size,
        )
        est[positions] = group.est
        delete[:, positions] = group.delete
    pseudo = n_blocks * est[None, :] - (n_blocks - 1) * delete
    coef = est / n_bar
    coef_se = np.sqrt(np.var(pseudo, axis=0, ddof=1) / n_blocks) / n_bar
    return BatchedSLDSCResult(
        labels=tuple(m.label for m in models),
        coef=coef,
        coef_se=coef_se,
        p=norm.sf(coef / coef_se),
    )


def _group_by_controls(
    models: Sequence[SLDSCModel],
) -> dict[tuple[int, ...], list[int]]:
    """Positions of the models in each control group, keyed by the group's control columns."""
    groups: dict[tuple[int, ...], list[int]] = defaultdict(list)
    for i, model in enumerate(models):
        groups[model.control_indices].append(i)
    return dict(groups)


def _shared_design(
    inputs: SLDSCRegressionInputs, controls: tuple[int, ...], scale: np.ndarray
) -> _SharedDesign:
    control_ld = inputs.annotation_ld[:, list(controls)].astype(np.float64)
    unscaled = np.hstack([control_ld, inputs.baseline_ld])
    return _SharedDesign(
        ld_sum=unscaled.sum(axis=1),
        m_sum=float(
            inputs.annotation_m[list(controls)].sum() + inputs.baseline_m.sum()
        ),
        design=np.hstack([unscaled * scale[:, None], np.ones((unscaled.shape[0], 1))]),
    )


def _exact_group(
    inputs: SLDSCRegressionInputs,
    shared: _SharedDesign,
    targets: list[int],
    scale: np.ndarray,
    separators: np.ndarray,
    chunk_size: int,
) -> _TargetEstimates:
    """Exact-weight target estimates of one control group, chunk_size models at a time."""
    upper_triangle = np.triu_indices(shared.design.shape[1])
    ests, deletes = [], []
    for start in range(0, len(targets), chunk_size):
        chunk = targets[start : start + chunk_size]
        unscaled_targets = inputs.annotation_ld[:, chunk].astype(np.float64)
        weights = ldsc_hsq_weights(
            ld_tot=shared.ld_sum[:, None] + unscaled_targets,
            w_ld=inputs.w_ld,
            n=inputs.n,
            m_tot=shared.m_sum + inputs.annotation_m[chunk],
            chi2=inputs.chi2,
        )
        blocks = _exact_block_normal_equations(
            shared_design=shared.design,
            weights=weights,
            scaled_targets=unscaled_targets * scale[:, None],
            chi2=inputs.chi2,
            separators=separators,
            upper_triangle=upper_triangle,
        )
        gram_tot = blocks.gram.sum(axis=0)
        rhs_tot = blocks.rhs.sum(axis=0)
        ests.append(_solve_first(gram=gram_tot, rhs=rhs_tot))
        deletes.append(
            _solve_first(
                gram=gram_tot[None] - blocks.gram, rhs=rhs_tot[None] - blocks.rhs
            )
        )
    return _TargetEstimates(
        est=np.concatenate(ests), delete=np.concatenate(deletes, axis=1)
    )


def _exact_block_normal_equations(
    shared_design: np.ndarray,
    weights: np.ndarray,
    scaled_targets: np.ndarray,
    chi2: np.ndarray,
    separators: np.ndarray,
    upper_triangle: tuple[np.ndarray, np.ndarray],
) -> _BlockNormalEquations:
    """
    Per-block weighted least-squares normal equations for a chunk of k models of one control
    group, each with its own weights.

    Parameters (S regression SNPs, q shared columns, k models; see the module docstring):
        shared_design: (S, q) matrix A of scaled control and baseline LD scores plus the
            intercept column, common to all k models.
        weights: (S, k); column c is the weight vector w_c of model c (W_c = diag(w_c)).
        scaled_targets: (S, k); column c is the scaled target LD score t_c of model c.
        chi2: (S,) response y.
        separators: (n_blocks + 1,) jackknife block boundaries; block b is rows
            separators[b]:separators[b + 1].
        upper_triangle: row and column indices of the upper triangle of a q x q matrix
            (np.triu_indices(q)).

    Model c has design X_c = [t_c, A] (S x p, p = q + 1) and weighted least-squares estimate
    beta_c solving the normal equations (X_c' W_c X_c) beta_c = X_c' W_c y. In blocks, with the
    target coefficient first:

        [ t_c' W_c t_c    t_c' W_c A ] [ beta_t ]   [ t_c' W_c y ]
        [ A' W_c t_c      A' W_c A   ] [ beta_A ] = [ A' W_c y   ]

    Every entry is a sum over SNPs, so it splits into per-block sums. For each block b this
    returns gram[b, c] = X_c' W_c X_c and rhs[b, c] = X_c' W_c y restricted to block b's SNPs,
    computed for all k models together:
        A' W_c A for every c: one product of the (s_b, q(q+1)/2) matrix of pairwise column
            products of A (upper triangle only) with the (s_b, k) weights, then mirrored.
        A' W_c t_c for every c: A' (weights * scaled_targets), one (q, s_b) x (s_b, k) product.
        t_c' W_c t_c, t_c' W_c y: column sums of weights * scaled_targets * scaled_targets and
            weights * scaled_targets * y.
        A' W_c y for every c: A' (weights * y), one (q, s_b) x (s_b, k) product.
    """
    q = shared_design.shape[1]
    k = weights.shape[1]
    n_blocks = len(separators) - 1
    gram = np.empty((n_blocks, k, q + 1, q + 1))
    rhs = np.empty((n_blocks, k, q + 1))
    diagonal = np.arange(q)
    rows_i, rows_j = upper_triangle
    for b in range(n_blocks):
        block = slice(separators[b], separators[b + 1])
        a_b = shared_design[block]
        w_b = weights[block]
        t_b = scaled_targets[block]
        y_b = chi2[block]
        weighted_targets = w_b * t_b
        upper = np.zeros((k, q, q))
        upper[:, rows_i, rows_j] = ((a_b[:, rows_i] * a_b[:, rows_j]).T @ w_b).T
        a_w_a = upper + np.transpose(upper, (0, 2, 1))
        a_w_a[:, diagonal, diagonal] /= 2.0
        a_w_t = (a_b.T @ weighted_targets).T  # (k, q)
        gram[b, :, 1:, 1:] = a_w_a
        gram[b, :, 0, 1:] = a_w_t
        gram[b, :, 1:, 0] = a_w_t
        gram[b, :, 0, 0] = (weighted_targets * t_b).sum(axis=0)
        rhs[b, :, 0] = (weighted_targets * y_b[:, None]).sum(axis=0)
        rhs[b, :, 1:] = (a_b.T @ (w_b * y_b[:, None])).T
    return _BlockNormalEquations(gram=gram, rhs=rhs)


def _solve_first(gram: np.ndarray, rhs: np.ndarray) -> np.ndarray:
    """
    First (target) coefficient of batched linear solves gram @ beta = rhs.
    gram: (..., p, p); rhs: (..., p); returns (...).
    """
    return np.linalg.solve(gram, rhs[..., None])[..., 0, 0]


def _shared_controls_group(
    inputs: SLDSCRegressionInputs,
    shared: _SharedDesign,
    targets: list[int],
    scale: np.ndarray,
    separators: np.ndarray,
    chunk_size: int,
) -> _TargetEstimates:
    """
    Shared-weight target estimates of one control group. Every model uses the weight w of the
    model without a target column, so W = diag(w) is common to the group: the per-block
    A' W A (q x q) and A' W y (q,) are formed once, the full-data and delete-block A' W A are
    factored once, and each chunk of k targets only adds
    A' W t_c, t_c' W t_c and t_c' W y (dimensions as in the module docstring).
    """
    w = ldsc_hsq_weights(
        ld_tot=shared.ld_sum[:, None],
        w_ld=inputs.w_ld,
        n=inputs.n,
        m_tot=np.array([shared.m_sum]),
        chi2=inputs.chi2,
    )[:, 0]
    a = shared.design  # (S, q)
    y = inputs.chi2  # (S)
    blocks = [
        slice(separators[b], separators[b + 1]) for b in range(len(separators) - 1)
    ]
    a_w_a = np.stack([a[r].T @ (w[r, None] * a[r]) for r in blocks])  # (n_blocks, q, q)
    a_w_y = np.stack([a[r].T @ (w[r] * y[r]) for r in blocks])  # (n_blocks, q)
    full_gram = _factor_shared_gram(a_w_a=a_w_a.sum(axis=0), a_w_y=a_w_y.sum(axis=0))
    delete_gram = _factor_shared_gram(
        a_w_a=a_w_a.sum(axis=0)[None] - a_w_a, a_w_y=a_w_y.sum(axis=0)[None] - a_w_y
    )
    ests, deletes = [], []
    for start in range(0, len(targets), chunk_size):
        chunk = targets[start : start + chunk_size]
        t = inputs.annotation_ld[:, chunk].astype(np.float64) * scale[:, None]  # (S, k)
        w_t = w[:, None] * t  # (S, k)
        a_w_t = np.stack([a[r].T @ w_t[r] for r in blocks])  # (n_blocks, q, k)
        t_w_t = np.stack([(w_t[r] * t[r]).sum(axis=0) for r in blocks])  # (n_blocks, k)
        t_w_y = np.stack(
            [(w_t[r] * y[r, None]).sum(axis=0) for r in blocks]
        )  # (n_blocks, k)
        ests.append(
            _fwl_target_coef(
                gram=full_gram,
                a_w_t=a_w_t.sum(axis=0),
                t_w_t=t_w_t.sum(axis=0),
                t_w_y=t_w_y.sum(axis=0),
            )
        )
        deletes.append(
            _fwl_target_coef(
                gram=delete_gram,
                a_w_t=a_w_t.sum(axis=0)[None] - a_w_t,
                t_w_t=t_w_t.sum(axis=0)[None] - t_w_t,
                t_w_y=t_w_y.sum(axis=0)[None] - t_w_y,
            )
        )
    return _TargetEstimates(
        est=np.concatenate(ests), delete=np.concatenate(deletes, axis=1)
    )


@frozen(slots=True)
class _FactoredSharedGram:
    """
    The shared Gram matrix A' W A of a control group, LU-factored once and reused by every chunk
    of targets, together with the coefficients of y regressed on A, which do not depend on the
    target either.

    Leading axes as in _fwl_target_coef (none for the full data, n_blocks for delete values):
        lu, piv: LU factorization of A' W A, as returned by scipy.linalg.lu_factor:
            lu (..., q, q), piv (..., q).
        y_coef: (..., q) (A' W A)^-1 A' W y.
    """

    lu: np.ndarray
    piv: np.ndarray
    y_coef: np.ndarray

    def __attrs_post_init__(self) -> None:
        assert self.lu.ndim >= 2 and self.lu.shape[-1] == self.lu.shape[-2], (
            f"lu must be (..., q, q), got {self.lu.shape}"
        )
        assert self.piv.shape == self.lu.shape[:-1], (
            f"piv shape {self.piv.shape} does not match lu shape {self.lu.shape}"
        )
        assert self.y_coef.shape == self.lu.shape[:-1], (
            f"y_coef shape {self.y_coef.shape} does not match lu shape {self.lu.shape}"
        )


def _factor_shared_gram(a_w_a: np.ndarray, a_w_y: np.ndarray) -> _FactoredSharedGram:
    """
    a_w_a: (..., q, q) A' W A; a_w_y: (..., q) A' W y.
    """
    lu, piv = lu_factor(a_w_a)
    return _FactoredSharedGram(
        lu=lu, piv=piv, y_coef=lu_solve((lu, piv), a_w_y[..., None])[..., 0]
    )


def _fwl_target_coef(
    gram: _FactoredSharedGram,
    a_w_t: np.ndarray,
    t_w_t: np.ndarray,
    t_w_y: np.ndarray,
) -> np.ndarray:
    """
    Target coefficient of the regressions y ~ [t_c, A] with a common weight matrix W, by the
    Frisch-Waugh-Lovell identity: residualize t_c and y on A, then regress residual on residual,

        beta_t = (t_c' W y - t_c' W A (A' W A)^-1 A' W y) / (t_c' W t_c - t_c' W A (A' W A)^-1 A' W t_c).

    Broadcasts over leading axes (none for the full-data estimate, n_blocks for delete values):
        gram: factored A' W A and (A' W A)^-1 A' W y.
        a_w_t: (..., q, k) A' W t_c for each of k targets.
        t_w_t, t_w_y: (..., k) t_c' W t_c and t_c' W y.
    Returns (..., k).
    """
    t_coef = lu_solve((gram.lu, gram.piv), a_w_t)  # (..., q, k)
    numerator = t_w_y - (a_w_t * gram.y_coef[..., None]).sum(axis=-2)
    denominator = t_w_t - (a_w_t * t_coef).sum(axis=-2)
    return numerator / denominator
