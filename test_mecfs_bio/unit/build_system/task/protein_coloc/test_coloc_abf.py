import math

import numpy as np
import pytest

from mecfs_bio.build_system.task.protein_coloc.coloc_abf import (
    CaseControlTraitModel,
    ColocAbf,
    ColocInput,
    ColocPriors,
    QuantitativeTraitModel,
    estimate_sd_y,
)

_N_VARIANTS = 200
_SE = 0.05


def _input(z1_at: dict[int, float], z2_at: dict[int, float]) -> ColocInput:
    z1 = np.zeros(_N_VARIANTS)
    z2 = np.zeros(_N_VARIANTS)
    for i, z in z1_at.items():
        z1[i] = z
    for i, z in z2_at.items():
        z2[i] = z
    se = np.full(_N_VARIANTS, _SE)
    return ColocInput(
        beta1=z1 * se,
        se1=se,
        beta2=z2 * se,
        se2=se.copy(),
        maf1=None,
        n1=None,
        maf2=None,
        n2=None,
        model1=QuantitativeTraitModel(sd_y=1.0),
        model2=QuantitativeTraitModel(sd_y=1.0),
    )


def test_shared_causal_variant_gives_high_h4() -> None:
    post = ColocAbf().run(_input(z1_at={50: 10.0}, z2_at={50: 10.0}))
    assert post.pp_h4 > 0.9


def test_distinct_causal_variants_give_high_h3() -> None:
    post = ColocAbf().run(_input(z1_at={50: 10.0}, z2_at={150: 10.0}))
    assert post.pp_h3 > 0.9


def test_no_signal_gives_high_h0() -> None:
    assert ColocAbf().run(_input(z1_at={}, z2_at={})).pp_h0 > 0.9


def test_posteriors_sum_to_one() -> None:
    p = ColocAbf().run(_input(z1_at={50: 6.0}, z2_at={52: 4.0}))
    assert abs(p.pp_h0 + p.pp_h1 + p.pp_h2 + p.pp_h3 + p.pp_h4 - 1) < 1e-12


def test_larger_p12_increases_h4() -> None:
    data = _input(z1_at={50: 5.0}, z2_at={50: 5.0})
    assert (
        ColocAbf(ColocPriors(p12=1e-4)).run(data).pp_h4
        > ColocAbf(ColocPriors(p12=1e-6)).run(data).pp_h4
    )


def _labf(z: float, se: float, w: float) -> float:
    r = w**2 / (w**2 + se**2)
    return 0.5 * (math.log(1 - r) + r * z**2)


def test_single_variant_matches_hand_computed_abf() -> None:
    se1, se2 = 0.04, 0.03
    z1, z2 = 4.0, 3.0
    data = ColocInput(
        beta1=np.array([z1 * se1]),
        se1=np.array([se1]),
        beta2=np.array([z2 * se2]),
        se2=np.array([se2]),
        maf1=None,
        n1=None,
        maf2=None,
        n2=None,
        model1=CaseControlTraitModel(case_fraction=0.3),
        model2=QuantitativeTraitModel(sd_y=1.0),
    )
    priors = ColocPriors()
    l1 = _labf(z1, se1, 0.2)
    l2 = _labf(z2, se2, 0.15)
    # With one variant H3 (two distinct causal variants) is impossible.
    weights = [
        1.0,
        priors.p1 * math.exp(l1),
        priors.p2 * math.exp(l2),
        priors.p12 * math.exp(l1 + l2),
    ]
    post = ColocAbf(priors).run(data)
    assert post.pp_h4 == pytest.approx(weights[3] / sum(weights), abs=1e-12)
    assert post.pp_h3 == 0.0


def test_sd_y_estimated_when_absent() -> None:
    rng = np.random.default_rng(0)
    n_variants = 2000
    sd_y = 2.0
    maf = rng.uniform(0.05, 0.5, n_variants)
    n = np.full(n_variants, 20_000.0)
    se = sd_y / np.sqrt(2 * n * maf * (1 - maf)) * rng.uniform(0.98, 1.02, n_variants)
    assert estimate_sd_y(se=se, maf=maf, n=n) == pytest.approx(sd_y, rel=0.05)


def test_quantitative_model_without_sd_y_requires_maf_and_n() -> None:
    data = _input(z1_at={}, z2_at={})
    with pytest.raises(AssertionError):
        ColocInput(
            beta1=data.beta1,
            se1=data.se1,
            beta2=data.beta2,
            se2=data.se2,
            maf1=None,
            n1=None,
            maf2=None,
            n2=None,
            model1=QuantitativeTraitModel(sd_y=None),
            model2=QuantitativeTraitModel(sd_y=1.0),
        )
