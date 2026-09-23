import numpy as np
import pytest

from mecfs_bio.build_system.task.annotation_weights.chromosome_blocked_ridge import (
    ChromRidgeBlock,
    accumulate_block,
    combine,
    fit,
    select_alpha_loco,
)


def _chrom_blocks(
    seed: int, beta: np.ndarray, n_chrom: int = 6, n_per: int = 200
) -> dict[int, ChromRidgeBlock]:
    rng = np.random.default_rng(seed)
    blocks = {}
    for chrom in range(1, n_chrom + 1):
        x = rng.normal(size=(n_per, beta.size))
        y = x @ beta + rng.normal(scale=0.01, size=n_per)
        blocks[chrom] = accumulate_block(x=x, y=y)
    return blocks


def test_fit_recovers_known_beta_unweighted():
    beta = np.array([2.0, -1.0, 0.5])
    blocks = _chrom_blocks(0, beta)
    result = fit(combine(list(blocks.values())), alpha=1e-6)
    assert np.allclose(result.beta_raw, beta, atol=1e-2)


def test_combine_is_additive():
    beta = np.array([1.0, 0.0])
    blocks = list(_chrom_blocks(1, beta, n_chrom=2).values())
    merged = combine(blocks)
    assert merged.sw == pytest.approx(blocks[0].sw + blocks[1].sw)
    assert np.allclose(merged.swxx, blocks[0].swxx + blocks[1].swxx)


def test_select_alpha_loco_prefers_small_alpha_on_clean_signal():
    beta = np.array([3.0, -2.0, 1.0])
    blocks = _chrom_blocks(2, beta)
    sel = select_alpha_loco(blocks, alphas=(1e-6, 1.0, 1e3, 1e6))
    assert sel.alpha == 1e-6
    assert sel.mean_r2 > 0.999


def test_weights_change_the_fit():
    # Half the points are corrupted; down-weighting them should recover beta
    # better than the unweighted fit does.
    rng = np.random.default_rng(3)
    x = rng.normal(size=(500, 2))
    beta = np.array([1.0, -1.0])
    y = x @ beta + rng.normal(scale=0.1, size=500)
    y[:250] += 5.0 * x[:250, 0]
    w = np.ones(500)
    w[:250] = 1e-6
    weighted = fit(accumulate_block(x=x, y=y, w=w), alpha=1e-6)
    unweighted = fit(accumulate_block(x=x, y=y), alpha=1e-6)
    assert np.linalg.norm(weighted.beta_raw - beta) < np.linalg.norm(
        unweighted.beta_raw - beta
    )
