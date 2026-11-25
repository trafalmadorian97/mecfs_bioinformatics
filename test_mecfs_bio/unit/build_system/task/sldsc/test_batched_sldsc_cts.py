"""
Tests for the batched S-LDSC kernel. The exact-mode oracle is ldsc's Hsq class
(Hsq = "h-squared", ldsc's heritability estimator in gwaslab/extension/ldsc/ldsc_regressions.py),
called the way cell_type_specific calls it.
"""

import numpy as np
import pytest
from attrs import frozen
from gwaslab.extension.ldsc.ldsc_regressions import Hsq

from mecfs_bio.build_system.task.sldsc.batched_sldsc_cts import (
    SLDSCModel,
    SLDSCRegressionInputs,
    batched_sldsc_cts,
    ldsc_block_separators,
    ldsc_hsq_weights,
)


@frozen(slots=True)
class _OracleFit:
    coef: float
    se: float


def _random_inputs(n_snps: int, seed: int = 0) -> SLDSCRegressionInputs:
    rng = np.random.default_rng(seed)
    baseline_ld = np.column_stack(
        [
            rng.uniform(1, 50, n_snps),
            rng.uniform(0, 10, n_snps),
            rng.uniform(0, 5, n_snps),
        ]
    )
    annotation_ld = rng.uniform(0, 10, (n_snps, 5)).astype(np.float32)
    n = np.full(n_snps, 40_000.0) + rng.integers(0, 2_000, n_snps)
    chi2 = rng.chisquare(1, n_snps) * (
        1.0 + 0.02 * baseline_ld[:, 0] + 0.1 * annotation_ld[:, 0]
    )
    return SLDSCRegressionInputs(
        chi2=chi2,
        n=n,
        w_ld=rng.uniform(1, 40, n_snps),
        baseline_ld=baseline_ld,
        baseline_m=np.array([5000.0, 1000.0, 500.0]),
        annotation_ld=annotation_ld,
        annotation_m=np.array([300.0, 400.0, 500.0, 2000.0, 1500.0]),
    )


# Two control groups: columns 3 and 4 are controls; 0, 1, 2 are targets.
_MODELS = (
    SLDSCModel("a", 0, (3,)),
    SLDSCModel("b", 1, (3,)),
    SLDSCModel("c", 2, (4,)),
    SLDSCModel("d", 0, (3, 4)),
)


def _gwaslab_oracle(
    inputs: SLDSCRegressionInputs, model: SLDSCModel, n_blocks: int
) -> _OracleFit:
    cols = [model.target_index, *model.control_indices]
    ref_ld = np.column_stack(
        [inputs.annotation_ld[:, cols].astype(np.float64), inputs.baseline_ld]
    )
    m = np.concatenate([inputs.annotation_m[cols], inputs.baseline_m]).reshape(1, -1)
    s = inputs.chi2.shape[0]
    hsq = Hsq(
        y=inputs.chi2.reshape(s, 1),
        x=ref_ld,
        w=inputs.w_ld.reshape(s, 1),
        N=inputs.n.reshape(s, 1),
        M=m,
        n_blocks=min(s, n_blocks),
        intercept=None,
        twostep=None,
        old_weights=True,
    )
    return _OracleFit(coef=float(hsq.coef[0]), se=float(hsq.coef_se[0]))


@pytest.mark.parametrize(
    "n_snps", [3_000, 150]
)  # 150 < 200 exercises n_blocks = min(S, 200)
def test_exact_matches_gwaslab_hsq(n_snps: int) -> None:
    inputs = _random_inputs(n_snps)
    result = batched_sldsc_cts(inputs, _MODELS, weighting="exact", chunk_size=2)
    for i, model in enumerate(_MODELS):
        oracle = _gwaslab_oracle(inputs, model, n_blocks=200)
        np.testing.assert_allclose(result.coef[i], oracle.coef, rtol=1e-8)
        np.testing.assert_allclose(result.coef_se[i], oracle.se, rtol=1e-8)
    assert result.labels == tuple(m.label for m in _MODELS)


def test_exact_is_independent_of_chunk_size() -> None:
    inputs = _random_inputs(3_000)
    small = batched_sldsc_cts(inputs, _MODELS, weighting="exact", chunk_size=1)
    large = batched_sldsc_cts(inputs, _MODELS, weighting="exact", chunk_size=64)
    np.testing.assert_allclose(small.coef, large.coef, rtol=1e-12)
    np.testing.assert_allclose(small.coef_se, large.coef_se, rtol=1e-12)


def _dense_shared_weight_oracle(
    inputs: SLDSCRegressionInputs, model: SLDSCModel, n_blocks: int
) -> _OracleFit:
    """Direct weighted least squares + leave-one-block-out refits, with the shared weight."""
    s = inputs.chi2.shape[0]
    controls = inputs.annotation_ld[:, list(model.control_indices)].astype(np.float64)
    shared_unscaled = np.hstack([controls, inputs.baseline_ld])
    m_shared = (
        inputs.annotation_m[list(model.control_indices)].sum() + inputs.baseline_m.sum()
    )
    w = ldsc_hsq_weights(
        ld_tot=shared_unscaled.sum(axis=1)[:, None],
        w_ld=inputs.w_ld,
        n=inputs.n,
        m_tot=np.array([m_shared]),
        chi2=inputs.chi2,
    )[:, 0]
    n_bar = inputs.n.mean()
    target = inputs.annotation_ld[:, model.target_index].astype(np.float64)
    x = np.column_stack(
        [
            np.column_stack([target, shared_unscaled]) * (inputs.n / n_bar)[:, None],
            np.ones(s),
        ]
    )
    sw = np.sqrt(w)

    def fit(rows: np.ndarray) -> float:
        beta, *_ = np.linalg.lstsq(
            x[rows] * sw[rows, None], inputs.chi2[rows] * sw[rows], rcond=None
        )
        return float(beta[0])

    nb = min(s, n_blocks)
    seps = ldsc_block_separators(s, nb)
    full = fit(np.arange(s))
    delete = np.array([fit(np.r_[0 : seps[j], seps[j + 1] : s]) for j in range(nb)])
    pseudo = nb * full - (nb - 1) * delete
    return _OracleFit(
        coef=full / n_bar, se=float(np.sqrt(np.var(pseudo, ddof=1) / nb) / n_bar)
    )


def test_shared_controls_matches_dense_oracle() -> None:
    inputs = _random_inputs(3_000)
    result = batched_sldsc_cts(
        inputs, _MODELS, weighting="shared_controls", chunk_size=2
    )
    for i, model in enumerate(_MODELS):
        oracle = _dense_shared_weight_oracle(inputs, model, n_blocks=200)
        np.testing.assert_allclose(result.coef[i], oracle.coef, rtol=1e-8)
        np.testing.assert_allclose(result.coef_se[i], oracle.se, rtol=1e-8)


def test_shared_controls_close_to_exact_on_random_data() -> None:
    # Sanity check, not a requirement: shared weights should not move z-scores far.
    inputs = _random_inputs(3_000)
    exact = batched_sldsc_cts(inputs, _MODELS, weighting="exact")
    shared = batched_sldsc_cts(inputs, _MODELS, weighting="shared_controls")
    np.testing.assert_allclose(
        shared.coef / shared.coef_se, exact.coef / exact.coef_se, atol=0.2
    )
