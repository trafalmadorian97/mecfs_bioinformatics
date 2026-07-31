"""
Validate a batched, vectorized LDSC heritability kernel against the repo's exact
per-protein implementation (GenomicSEM port, genomic_sem_ldsc.estimate_h2).

Motivation (see plan): all PPP proteins share one variant-index row order, so the
index<->LD-score alignment, ld/wLD/M and the jackknife blocks are protein-invariant.
That lets us process K proteins at once: stack chi^2 into an (n_snps x K) matrix and
run the weighted LD-score regression + block jackknife for all K in vectorized numpy.
Per-protein filtering (missing variants + chi^2 filter) is encoded as ZERO regression
weights instead of dropping rows, so every protein keeps the same shared row set and
the same shared contiguous jackknife blocks.

Key expectations this probe checks:
  - h^2 POINT estimate is machine-precision identical (zero-weight == dropped, exactly,
    since blocks only partition the weighted normal-equation sums).
  - h^2 jackknife SE differs only marginally: the exact method cuts the *kept* SNPs into
    n_blocks contiguous blocks, whereas the batched method cuts the *full* shared set
    (filtered rows zero-weighted) into n_blocks blocks, so block membership shifts by a
    few SNPs. We quantify that SE discrepancy here.

Run:
    pixi r python experiments/claude/ppp_ldsc/batched_vs_exact_h2_probe.py

Reads only local asset-store files; writes a log to experiments/claude/ppp_ldsc/logs/.
"""

from __future__ import annotations

import glob
import os
import time
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_ldsc import (
    estimate_h2,
)

REPO = Path(__file__).resolve().parents[3]
STORE = REPO / "assets" / "base_asset_store"
INDEX_PATH = (
    STORE
    / "reference_data/ukbb_ppp_variant_index/hapmap_3_membership_list/processed/ppp_variant_index.parquet"
)
LD_DIR = (
    STORE
    / "reference_data/linkage_disequilibrium_scores/thousand_genomes_phase_3_v1/extracted/thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)
SLIM_GLOB = str(
    STORE / "gwas/ukbb_ppp/*/aligned/hapmap_3_index/ppp_slim_hapmap_3_*.parquet.zstd"
)
LOG_DIR = REPO / "experiments/claude/ppp_ldsc/logs"

N_BLOCKS = 200
N_PROTEINS = 15
N_CHROM = 22


def read_ld_scores(ld_dir: Path, n_chrom: int) -> tuple[pl.DataFrame, float]:
    """Read LDscore.<chr>.l2.ldscore.gz (SNP, CHR, BP, L2) and .l2.M_5_50, summing M."""
    frames = []
    m_total = 0.0
    for chrom in range(1, n_chrom + 1):
        score = pl.read_csv(
            ld_dir / f"LDscore.{chrom}.l2.ldscore.gz",
            separator="\t",
        ).select("CHR", "SNP", "BP", "L2")
        frames.append(score)
        m_val = pl.read_csv(
            ld_dir / f"LDscore.{chrom}.l2.M_5_50", separator="\t", has_header=False
        )
        m_total += float(m_val.to_numpy().sum())
    return pl.concat(frames), m_total


def build_context() -> dict:
    """Shared, protein-invariant context: the regression SNP set (index rows that carry
    an LD score), genome-sorted, with ld/wLD, index row positions, EAF, and M."""
    index = pl.read_parquet(INDEX_PATH).with_row_index("__row__")
    ld_df, m_total = read_ld_scores(LD_DIR, N_CHROM)

    # rsID -> L2 join (inner). Keep index (genomic) order, which is already sorted by
    # (CHR, POS, EA, NEA) as ConstructPppVariantIndexTask produces it.
    merged = (
        index.join(ld_df, left_on="rsID", right_on="SNP", how="inner")
        .sort(["CHR", "POS"])  # contiguous blocks need genome order
    )
    return {
        "row_pos": merged["__row__"].to_numpy(),  # positions into the slim files
        "ld": merged["L2"].to_numpy().astype(float),
        "eaf": merged["EAF"].to_numpy().astype(float),
        "m": m_total,
        "n_snps": merged.height,
    }


def effective_n(se: np.ndarray, eaf: np.ndarray) -> float:
    """Per-protein N ~= median 1/(2 MAF (1-MAF) SE^2) (INT trait, var~=1)."""
    maf = np.minimum(eaf, 1 - eaf)
    ok = np.isfinite(se) & (se > 0) & (maf > 0)
    return float(np.median(1.0 / (2 * maf[ok] * (1 - maf[ok]) * se[ok] ** 2)))


def load_protein_chi2(path: str, ctx: dict) -> tuple[np.ndarray, float]:
    """chi^2 = (BETA/SE)^2 at the shared context rows (NaN where the protein is missing),
    plus the protein's effective N."""
    df = pl.read_parquet(path)
    beta = df["BETA"].to_numpy().astype(float)[ctx["row_pos"]]
    se = df["SE"].to_numpy().astype(float)[ctx["row_pos"]]
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = (beta / se) ** 2
    n = effective_n(se, ctx["eaf"])
    return chi2, n


# --------------------------------------------------------------------------------------
# Exact reference (repo code) + its jackknife SE, mirroring run_ldsc's V computation.
# --------------------------------------------------------------------------------------
def exact_h2(chi2: np.ndarray, ld: np.ndarray, n_scalar: float, m: float) -> tuple[float, float]:
    keep = np.isfinite(chi2)
    thr = max(0.001 * n_scalar, 80.0)
    keep &= chi2 <= thr
    chi = chi2[keep]
    ldk = ld[keep]
    n = np.full(chi.shape, n_scalar)
    est = estimate_h2(chi=chi, ld_raw=ldk, wld_raw=ldk, n=n, m=m, n_blocks=N_BLOCKS)
    # run_ldsc: V_h2 = var(pseudo, ddof=1) / (n_bar*sqrt(n_blocks)/m)^2 ; se = sqrt(V).
    denom = (est.n_bar * np.sqrt(N_BLOCKS) / m) ** 2
    se = float(np.sqrt(np.var(est.pseudo_coef, ddof=1) / denom))
    return est.reg_tot, se


# --------------------------------------------------------------------------------------
# Batched vectorized kernel: K proteins on the shared row set, drops as zero weights,
# shared contiguous blocks over the FULL set.
# --------------------------------------------------------------------------------------
def block_bounds(n_snps: int, n_blocks: int) -> list[tuple[int, int]]:
    sf = np.floor(np.linspace(1, n_snps, n_blocks + 1)).astype(int)
    out = []
    for p in range(n_blocks):
        a = int(sf[p]) - 1
        b = int(sf[p + 1]) - 1 if p < n_blocks - 1 else n_snps
        out.append((a, b))
    return out


def batched_h2(
    chi2: np.ndarray, ld: np.ndarray, n: np.ndarray, m: float
) -> tuple[np.ndarray, np.ndarray]:
    """chi2: (S, K) with NaN for missing; ld: (S,); n: (K,) per-protein N; returns
    (h2[K], se[K])."""
    s, k = chi2.shape
    thr = np.maximum(0.001 * n, 80.0)  # (K,)
    keep = np.isfinite(chi2) & (chi2 <= thr[None, :])  # (S, K)
    chi = np.where(keep, chi2, 0.0)

    # tot.agg per protein over KEPT SNPs: m*(mean(chi)-1)/mean(ld*N).
    cnt = keep.sum(0).astype(float)  # (K,)
    mean_chi = chi.sum(0) / cnt
    ld_col = ld[:, None]
    mean_ldn = (np.where(keep, ld_col, 0.0).sum(0) * n) / cnt  # mean(ld)*N == mean(ld*N)
    tot_agg = np.clip(m * (mean_chi - 1.0) / mean_ldn, 0.0, 1.0)  # (K,)

    # WLS weight omega = het*oc, zeroed on non-kept. NOTE: estimate_h2 multiplies BOTH
    # design and response by initial_w = sqrt(het*oc), so the effective normal-equation
    # weight is (sqrt(het*oc))^2 = het*oc. The sum-to-1 normalization it applies is a
    # common scalar that cancels in solve(XtWX, XtWy), so we omit it here.
    ldm = np.maximum(ld, 1.0)[:, None]  # wLD == LD here
    c = (tot_agg * n / m)[None, :]  # (1,K)
    het = 1.0 / (2.0 * (1.0 + c * ldm) ** 2)
    oc = 1.0 / ldm
    w = het * oc * keep  # (S,K), effective WLS weight, zero where dropped

    # Per-block weighted sums for the 2x2 normal equations of design [ld, 1]:
    #   a=sum w*ld^2, b=sum w*ld, d=sum w ; e=sum w*ld*chi, f=sum w*chi.
    ld2 = (ld * ld)[:, None]
    bounds = block_bounds(s, N_BLOCKS)
    blk = np.empty((N_BLOCKS, k, 5))  # a,b,d,e,f
    for i, (lo, hi) in enumerate(bounds):
        wb = w[lo:hi]  # (rows,K)
        cb = chi[lo:hi]
        blk[i, :, 0] = (wb * ld2[lo:hi]).sum(0)
        blk[i, :, 1] = (wb * ld_col[lo:hi]).sum(0)
        blk[i, :, 2] = wb.sum(0)
        blk[i, :, 3] = (wb * ld_col[lo:hi] * cb).sum(0)
        blk[i, :, 4] = (wb * cb).sum(0)

    tot = blk.sum(0)  # (K,5)

    def slope(a, b, d, e, f):  # closed-form 2x2 solve, slope component only
        det = a * d - b * b
        return (d * e - b * f) / det

    slope_full = slope(tot[:, 0], tot[:, 1], tot[:, 2], tot[:, 3], tot[:, 4])  # (K,)

    # Leave-one-block-out slopes -> pseudo-values of the raw slope.
    loo = tot[None, :, :] - blk  # (n_blocks, K, 5)
    slope_loo = slope(
        loo[:, :, 0], loo[:, :, 1], loo[:, :, 2], loo[:, :, 3], loo[:, :, 4]
    )  # (n_blocks, K)
    pseudo = N_BLOCKS * slope_full[None, :] - (N_BLOCKS - 1) * slope_loo  # (n_blocks,K)

    h2 = slope_full / n * m
    denom = (n * np.sqrt(N_BLOCKS) / m) ** 2
    se = np.sqrt(np.var(pseudo, axis=0, ddof=1) / denom)
    return h2, se


def main() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = open(LOG_DIR / "batched_vs_exact_h2_probe.log", "w")

    def emit(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        log.write(line + "\n")

    t0 = time.time()
    ctx = build_context()
    emit(f"context: n_snps={ctx['n_snps']} M={ctx['m']:.0f} build={time.time()-t0:.1f}s")

    files = sorted(glob.glob(SLIM_GLOB))[:N_PROTEINS]
    chi2_cols, n_list, names = [], [], []
    for f in files:
        chi2, n = load_protein_chi2(f, ctx)
        chi2_cols.append(chi2)
        n_list.append(n)
        names.append(
            os.path.basename(f).replace("ppp_slim_hapmap_3_", "").replace(
                ".parquet.zstd", ""
            )
        )
    chi2_mat = np.column_stack(chi2_cols)  # (S, K)
    n_arr = np.array(n_list)

    # Exact (per-protein loop, repo estimate_h2).
    t1 = time.time()
    ex_h2, ex_se = [], []
    for j in range(len(files)):
        h2, se = exact_h2(chi2_mat[:, j], ctx["ld"], n_arr[j], ctx["m"])
        ex_h2.append(h2)
        ex_se.append(se)
    t_exact = time.time() - t1
    ex_h2 = np.array(ex_h2)
    ex_se = np.array(ex_se)

    # Batched (all K at once).
    t2 = time.time()
    ba_h2, ba_se = batched_h2(chi2_mat, ctx["ld"], n_arr, ctx["m"])
    t_batched = time.time() - t2

    emit("")
    emit(
        f"{'protein':28s} {'N':>7s} {'h2_exact':>10s} {'h2_batch':>10s} "
        f"{'dh2':>10s} {'se_exact':>9s} {'se_batch':>9s} {'dse_rel':>8s}"
    )
    for j, nm in enumerate(names):
        dh2 = ba_h2[j] - ex_h2[j]
        dse_rel = (ba_se[j] - ex_se[j]) / ex_se[j] if ex_se[j] else float("nan")
        emit(
            f"{nm:28s} {n_arr[j]:7.0f} {ex_h2[j]:10.5f} {ba_h2[j]:10.5f} "
            f"{dh2:10.2e} {ex_se[j]:9.5f} {ba_se[j]:9.5f} {dse_rel:8.2%}"
        )

    emit("")
    emit(f"max |dh2| (point estimate)      : {np.max(np.abs(ba_h2 - ex_h2)):.3e}")
    emit(f"max relative |dse| (jackknife SE): {np.max(np.abs((ba_se-ex_se)/ex_se)):.3%}")
    emit(
        f"timing: exact loop = {t_exact:.2f}s ({t_exact/len(files)*1000:.0f} ms/protein); "
        f"batched = {t_batched:.2f}s ({t_batched/len(files)*1000:.0f} ms/protein); "
        f"speedup = {t_exact/t_batched:.1f}x"
    )
    log.close()


if __name__ == "__main__":
    main()
