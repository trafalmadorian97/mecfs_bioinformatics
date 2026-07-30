"""
End-to-end validation of the batched cross-trait rg path against the GenomicSEM port.

This exercises the FULL trait side of PppProteinCrossTraitRgTask -- align_trait_to_context (with
deliberately strand-swapped trait rows, so the allele flip is on the hot path) feeding
estimate_trait_context + batched_rg -- and compares the resulting rg (trait vs each protein) to
genomic_sem_ldsc.run_ldsc run on the same simulated data written out as munged .sumstats files and
a reference LD-score directory.

Because run_ldsc computes each diagonal h2 on its own kept set and each covariance on the pair
intersection -- exactly the convention the batched kernel follows -- the rg POINT estimates should
agree to ~machine precision (the shared-block jackknife only perturbs the SE). All munged files are
written on the index effect-allele orientation, so run_ldsc's own pairwise flip is a no-op and the
two paths share one reference orientation.

Run: pixi r python experiments/claude/ppp_ldsc/batched_rg_vs_run_ldsc_probe.py
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_rg import (
    batched_rg,
    estimate_trait_context,
)
from mecfs_bio.build_system.task.ppp_ldsc.ppp_ldsc_context import build_ppp_ldsc_context
from mecfs_bio.build_system.task.ppp_ldsc.trait_alignment import align_trait_to_context
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_ldsc import run_ldsc
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ppp_database_constants import PPP_INDEX_IS_STRAND_AMBIGUOUS_COL

_LOG = Path(__file__).resolve().parent / "logs" / "batched_rg_vs_run_ldsc_probe.log"

N_SNPS = 6000
K = 3
N_SAMPLE = 20_000.0
N_BLOCKS = 200
SEED = 0


def _simulate() -> tuple[pl.DataFrame, float, np.ndarray, np.ndarray, np.ndarray]:
    """Return (index_df, m_total, ld, z_trait, z_protein) on the index effect-allele orientation."""
    rng = np.random.default_rng(SEED)
    ld = rng.uniform(1.0, 40.0, size=N_SNPS)
    shared = rng.standard_normal(N_SNPS) * np.sqrt(0.03 * ld)
    z_trait = shared + rng.standard_normal(N_SNPS) * np.sqrt(0.03 * ld) + rng.standard_normal(N_SNPS)
    z_protein = np.empty((N_SNPS, K))
    for j in range(K):
        z_protein[:, j] = (
            (0.4 + 0.25 * j) * shared
            + rng.standard_normal(N_SNPS) * np.sqrt(0.03 * ld)
            + rng.standard_normal(N_SNPS)
        )
    index_df = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: np.ones(N_SNPS, dtype=np.int64),
            GWASLAB_POS_COL: np.arange(1, N_SNPS + 1, dtype=np.int64),
            GWASLAB_RSID_COL: [f"rs{i}" for i in range(N_SNPS)],
            GWASLAB_EFFECT_ALLELE_COL: ["A"] * N_SNPS,
            GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * N_SNPS,
            PPP_INDEX_IS_STRAND_AMBIGUOUS_COL: [False] * N_SNPS,
        }
    )
    m_total = float(N_SNPS)
    return index_df, m_total, ld, z_trait, z_protein


def _write_munged(path: Path, rsids: list[str], z: np.ndarray) -> None:
    frame = pl.DataFrame(
        {"SNP": rsids, "N": np.full(z.shape[0], N_SAMPLE), "Z": z, "A1": ["A"] * z.shape[0], "A2": ["G"] * z.shape[0]}
    )
    frame.write_csv(path, separator="\t")


def _write_ld_dir(ld_dir: Path, rsids: list[str], ld: np.ndarray) -> None:
    ld_dir.mkdir(parents=True, exist_ok=True)
    frame = pl.DataFrame(
        {"CHR": np.ones(len(rsids), dtype=np.int64), "SNP": rsids, "BP": np.arange(1, len(rsids) + 1), "L2": ld}
    )
    with gzip.open(ld_dir / "1.l2.ldscore.gz", "wb") as handle:
        handle.write(frame.write_csv(separator="\t").encode())
    (ld_dir / "1.l2.M_5_50").write_text(f"{N_SNPS}\n")


def _run_ldsc_rg(tmp: Path, rsids: list[str], z_trait: np.ndarray, z_protein: np.ndarray) -> np.ndarray:
    munged = []
    trait_path = tmp / "trait.sumstats"
    _write_munged(trait_path, rsids, z_trait)
    munged.append(trait_path)
    for j in range(K):
        p = tmp / f"protein_{j}.sumstats"
        _write_munged(p, rsids, z_protein[:, j])
        munged.append(p)
    ld_dir = tmp / "ld"
    _write_ld_dir(ld_dir, rsids, _LD_CACHE)
    result = run_ldsc(
        munged_paths=munged,
        ld_dir=ld_dir,
        sample_prev=[None] * (K + 1),
        population_prev=[None] * (K + 1),
        n_blocks=N_BLOCKS,
        n_chrom=1,
    )
    assert result.S_Stand is not None
    return np.array([result.S_Stand[0, j + 1] for j in range(K)])


def _our_rg(index_df: pl.DataFrame, m_total: float, ld: np.ndarray, z_trait: np.ndarray, z_protein: np.ndarray) -> np.ndarray:
    context = build_ppp_ldsc_context(
        index_df,
        pl.DataFrame({"SNP": index_df[GWASLAB_RSID_COL], "L2": ld}),
        m_total,
        drop_strand_ambiguous=True,
        exclude_mhc=False,
    )
    context_variants = index_df.select(
        GWASLAB_RSID_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL
    ).gather(pl.Series(context.row_pos))

    # Build a trait dataframe with HALF the rows strand-swapped (EA/NEA flipped + BETA negated),
    # so alignment must flip them back; SE = 1 so z = BETA/SE = the canonical trait z.
    rng = np.random.default_rng(1)
    swap = rng.random(N_SNPS) < 0.5
    ea = np.where(swap, "G", "A")
    nea = np.where(swap, "A", "G")
    beta = np.where(swap, -z_trait, z_trait)
    trait_df = pl.DataFrame(
        {
            GWASLAB_RSID_COL: index_df[GWASLAB_RSID_COL],
            GWASLAB_EFFECT_ALLELE_COL: ea,
            GWASLAB_NON_EFFECT_ALLELE_COL: nea,
            GWASLAB_BETA_COL: beta,
            GWASLAB_SE_COL: np.ones(N_SNPS),
            GWASLAB_SAMPLE_SIZE_COLUMN: np.full(N_SNPS, N_SAMPLE),
        }
    )
    aligned = align_trait_to_context(
        trait_df, context_variants, trait_total_sample_size=None, min_trait_snps=100
    )
    # The alignment must recover the canonical trait z at the context SNPs.
    canonical = z_trait[context.row_pos]
    assert np.allclose(aligned.z, canonical), "alignment did not recover the canonical trait z"

    trait_ctx = estimate_trait_context(aligned.z, aligned.n, context.ld, context.m, n_blocks=N_BLOCKS)
    res = batched_rg(
        trait_ctx, z_protein[context.row_pos], context.ld, np.full(K, N_SAMPLE), context.m, n_blocks=N_BLOCKS
    )
    return res.rg


_LD_CACHE: np.ndarray  # set in main


def main() -> None:
    global _LD_CACHE
    index_df, m_total, ld, z_trait, z_protein = _simulate()
    _LD_CACHE = ld
    rsids = index_df[GWASLAB_RSID_COL].to_list()

    import tempfile

    with tempfile.TemporaryDirectory() as td:
        rg_run = _run_ldsc_rg(Path(td), rsids, z_trait, z_protein)
    rg_ours = _our_rg(index_df, m_total, ld, z_trait, z_protein)

    lines = ["protein  run_ldsc_rg   batched_rg    abs_diff"]
    for j in range(K):
        lines.append(f"{j:^7d}  {rg_run[j]:+.8f}  {rg_ours[j]:+.8f}  {abs(rg_run[j] - rg_ours[j]):.2e}")
    max_diff = float(np.max(np.abs(rg_run - rg_ours)))
    lines.append(f"\nmax |rg_run - rg_ours| = {max_diff:.2e}")
    lines.append("PASS" if max_diff < 1e-6 else "FAIL")
    report = "\n".join(lines)
    print(report)
    _LOG.write_text(report + "\n")
    assert max_diff < 1e-6, f"rg mismatch vs run_ldsc: {max_diff}"


if __name__ == "__main__":
    sys.exit(main())
