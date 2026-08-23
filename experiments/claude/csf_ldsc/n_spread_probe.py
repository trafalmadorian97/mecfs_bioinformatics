"""
Decide, on data, whether the CSF pQTL heritability task should use one scalar sample
size per aptamer (the median N) or the genuine per-SNP N.

Unlike UKB-PPP (N constant per protein), the Western et al. 2024 CSF sumstats carry a
per-variant N, so within a single aptamer N varies across its HapMap3 SNPs. Two things
decide whether that variation matters:

  1. How wide is the within-aptamer N distribution over the SNPs the LDSC regression
     actually uses (the LD-score context, strand-ambiguous SNPs dropped)?

  2. How much does the h2 point estimate actually move if we feed the estimator the real
     per-SNP N instead of the median?

Question 2 is the one that matters, and we can answer it exactly with the repo's already
validated reference kernel (genomic_sem_ldsc.estimate_h2), which accepts n as a per-SNP
array. Note how that estimator uses n: the design column stays ld (NOT ld*n), n enters
only the heteroscedasticity weights and the aggregate, and the final scaling is
reg_tot = slope / mean(n) * m. So "per-SNP N" vs "median N" here is dominantly
mean(N) vs median(N) in that scaling -- h2_persnp / h2_median ~ median(N) / mean(N) to
first order. We hold the kept-SNP set fixed (median-based chi^2 threshold) across both
runs so the ONLY thing that differs is the N handed to estimate_h2, isolating its effect.

If the within-aptamer N spread is tight and |h2 shift| is negligible, use the median and
reuse the PPP machinery unchanged. If it is wide, the per-SNP path (a batched-kernel
rewrite) is worth it.

Run:
    pixi r python experiments/claude/csf_ldsc/n_spread_probe.py

Reads only local asset-store files; writes a log to experiments/claude/csf_ldsc/logs/.
"""

from __future__ import annotations

import glob
import os
import random
import time
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_ldsc import (
    estimate_h2,
)
from mecfs_bio.constants.csf_database_constants import (
    CSF_INDEX_IS_STRAND_AMBIGUOUS_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

REPO = Path(__file__).resolve().parents[3]
STORE = REPO / "assets" / "base_asset_store"
INDEX_PATH = (
    STORE
    / "reference_data/csf_pqtl_variant_index/hapmap_3_membership_list/processed/csf_variant_index.parquet"
)
LD_DIR = (
    STORE
    / "reference_data/linkage_disequilibrium_scores/thousand_genomes_phase_3_v1/extracted/thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)
SLIM_GLOB = str(
    STORE / "gwas/western_csf_pqtl/*/aligned/hapmap_3_index/csf_slim_hapmap_3_*.parquet"
)
LOG_DIR = REPO / "experiments/claude/csf_ldsc/logs"

N_BLOCKS = 200
N_CHROM = 22
# A sample large enough to characterize the across-aptamer distribution without reading
# all 7,008 million-row files (~3 s each). Seeded for reproducibility.
N_APTAMERS = 150
SEED = 0


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
    """The regression SNP set: CSF index rows that carry an LD score, strand-ambiguous
    SNPs dropped (matching build_batched_ldsc_context's drop_strand_ambiguous=True), genome
    sorted for the block jackknife. Returns row positions into the slim files, ld and M."""
    index = pl.read_parquet(INDEX_PATH).with_row_index("__row__")
    index = index.filter(~pl.col(CSF_INDEX_IS_STRAND_AMBIGUOUS_COL))
    ld_df, m_total = read_ld_scores(LD_DIR, N_CHROM)
    merged = index.join(
        ld_df, left_on=GWASLAB_RSID_COL, right_on="SNP", how="inner"
    ).sort([GWASLAB_CHROM_COL, GWASLAB_POS_COL])
    return {
        "row_pos": merged["__row__"].to_numpy(),
        "ld": merged["L2"].to_numpy().astype(float),
        "m": m_total,
        "n_snps": merged.height,
    }


def load_protein(path: str, ctx: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """chi^2, per-SNP N and ld at the context rows (rows where the aptamer is present)."""
    df = pl.read_parquet(
        path, columns=[GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_SAMPLE_SIZE_COLUMN]
    )
    beta = df[GWASLAB_BETA_COL].to_numpy().astype(float)[ctx["row_pos"]]
    se = df[GWASLAB_SE_COL].to_numpy().astype(float)[ctx["row_pos"]]
    n = df[GWASLAB_SAMPLE_SIZE_COLUMN].to_numpy().astype(float)[ctx["row_pos"]]
    with np.errstate(invalid="ignore", divide="ignore"):
        chi2 = (beta / se) ** 2
    return chi2, n, ctx["ld"]


def h2_for_n(chi: np.ndarray, ld: np.ndarray, n: np.ndarray, m: float) -> float:
    """reg_tot from the repo's exact estimator over an already-filtered kept set."""
    est = estimate_h2(
        chi=chi, ld_raw=ld, wld_raw=ld, n=n, m=m, n_blocks=N_BLOCKS
    )
    return est.reg_tot


def summarize(values: np.ndarray) -> str:
    q = np.nanpercentile(values, [50, 75, 90, 95, 99, 100])
    return (
        f"median={q[0]:.4g} p75={q[1]:.4g} p90={q[2]:.4g} "
        f"p95={q[3]:.4g} p99={q[4]:.4g} max={q[5]:.4g}"
    )


def main() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = open(LOG_DIR / "n_spread_probe.log", "w")

    def emit(*a):
        line = " ".join(str(x) for x in a)
        print(line)
        log.write(line + "\n")

    t0 = time.time()
    ctx = build_context()
    emit(
        f"context: regression SNPs={ctx['n_snps']} M={ctx['m']:.0f} "
        f"build={time.time() - t0:.1f}s"
    )

    files = sorted(glob.glob(SLIM_GLOB))
    emit(f"slim files on disk: {len(files)}; sampling {N_APTAMERS} (seed={SEED})")
    random.seed(SEED)
    sample = sorted(random.sample(files, min(N_APTAMERS, len(files))))

    # Per-aptamer records.
    spread_iqr = []  # IQR(N)/median(N) over kept SNPs
    spread_p5_p95 = []  # (p95-p5)/median
    mean_over_median = []  # mean(N)/median(N): the first-order h2 ratio driver
    min_over_median = []
    max_over_median = []
    h2_rel_shift = []  # (h2_persnp - h2_median)/h2_median
    detail_rows = []

    t1 = time.time()
    for path in sample:
        chi2, n_all, ld_all = load_protein(path, ctx)
        median_n = float(np.nanmedian(n_all[np.isfinite(n_all)]))
        # Fixed kept set (median-based threshold) so only the N fed to the estimator differs.
        keep = np.isfinite(chi2) & np.isfinite(n_all)
        keep &= chi2 <= max(0.001 * median_n, 80.0)
        if keep.sum() < N_BLOCKS:
            continue
        chi = chi2[keep]
        ld = ld_all[keep]
        n_snp = n_all[keep]
        median_n = float(np.median(n_snp))

        q5, q25, q75, q95 = np.percentile(n_snp, [5, 25, 75, 95])
        spread_iqr.append((q75 - q25) / median_n)
        spread_p5_p95.append((q95 - q5) / median_n)
        mean_over_median.append(float(np.mean(n_snp)) / median_n)
        min_over_median.append(float(n_snp.min()) / median_n)
        max_over_median.append(float(n_snp.max()) / median_n)

        h2_median = h2_for_n(chi, ld, np.full(chi.shape, median_n), ctx["m"])
        h2_persnp = h2_for_n(chi, ld, n_snp, ctx["m"])
        rel = (h2_persnp - h2_median) / h2_median if h2_median != 0 else float("nan")
        h2_rel_shift.append(rel)

        detail_rows.append(
            (
                os.path.basename(path)
                .replace("csf_slim_hapmap_3_", "")
                .replace(".parquet", ""),
                keep.sum(),
                median_n,
                float(np.mean(n_snp)),
                (q75 - q25) / median_n,
                h2_median,
                h2_persnp,
                rel,
            )
        )

    emit(
        f"processed {len(detail_rows)} aptamers in {time.time() - t1:.1f}s "
        f"({(time.time() - t1) / max(len(detail_rows), 1):.2f}s each)"
    )

    # A few example rows (widest N spread first) so the log is legible.
    emit("")
    emit(
        f"{'aptamer':22s} {'kept':>7s} {'medN':>8s} {'meanN':>8s} "
        f"{'iqr/med':>8s} {'h2_med':>9s} {'h2_snp':>9s} {'h2_rel':>8s}"
    )
    for row in sorted(detail_rows, key=lambda r: -r[4])[:20]:
        nm, kept, medn, meann, iqr, h2m, h2s, rel = row
        emit(
            f"{nm:22s} {kept:7d} {medn:8.0f} {meann:8.0f} {iqr:8.2%} "
            f"{h2m:9.5f} {h2s:9.5f} {rel:8.2%}"
        )

    # Across-aptamer distributions -- the headline numbers for the decision.
    emit("")
    emit("=== within-aptamer N spread over the regression SNPs (across aptamers) ===")
    emit(f"IQR(N)/median(N)   : {summarize(np.array(spread_iqr))}")
    emit(f"(p95-p5)/median(N) : {summarize(np.array(spread_p5_p95))}")
    emit(f"mean(N)/median(N)  : {summarize(np.abs(np.array(mean_over_median) - 1.0))} "
         "(as |ratio-1|)")
    emit(f"min(N)/median(N)   : {summarize(np.array(min_over_median))}")
    emit(f"max(N)/median(N)   : {summarize(np.array(max_over_median))}")
    emit("")
    emit("=== |h2 shift| from median-N -> per-SNP-N (the decision-relevant number) ===")
    abs_shift = np.abs(np.array(h2_rel_shift))
    emit(f"|relative h2 shift|: {summarize(abs_shift)}")
    emit(
        f"aptamers with |h2 shift| > 1%: "
        f"{float(np.mean(abs_shift > 0.01)):.1%}; > 5%: "
        f"{float(np.mean(abs_shift > 0.05)):.1%}"
    )
    emit("")
    emit("INTERPRETATION:")
    emit("  If the |h2 shift| p95 is well under ~1-2%, the median is safe and we reuse the")
    emit("  PPP batched kernel unchanged. If it is large, the per-SNP path is worth the")
    emit("  batched-kernel rewrite.")
    log.close()


if __name__ == "__main__":
    main()
