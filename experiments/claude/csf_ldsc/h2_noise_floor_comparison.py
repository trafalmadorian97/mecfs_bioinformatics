"""
Compare the LDSC heritability diagnostics of the CSF pQTL database (Western et al. 2024,
median N ~ 3.3k) against UKB-PPP plasma (median N ~ 33.5k), to explain why the CSF table
has many large-magnitude NEGATIVE h2 estimates.

The headline: negative h2 is the LDSC noise floor, not a bug. LDSC h2 is an unconstrained
regression slope (chi2 on N*L/M); with little polygenic signal it scatters symmetrically
around zero. At CSF's ~10x smaller N the genome-wide inflation (mean_chi2) sits at the
null, so ~half the estimates go negative and none survive multiple-testing correction --
whereas PPP, with 10x the power, shows clearly positive h2 and >1,400 Bonferroni hits (all
positive). Negative h2 depends only on chi2 = (beta/se)^2, which is invariant to allele
orientation, so an alignment/sign error cannot produce it.

Run:
    pixi r python experiments/claude/csf_ldsc/h2_noise_floor_comparison.py

Reads only the two committed heritability parquet outputs; writes a log to
experiments/claude/csf_ldsc/logs/.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import polars as pl
import scipy.stats
from attrs import frozen

REPO = Path(__file__).resolve().parents[3]
STORE = REPO / "assets" / "base_asset_store"
CSF_PATH = (
    STORE / "gwas/western_csf_pqtl/csf_heritability/analysis/csf_heritability_hapmap_3.parquet"
)
PPP_PATH = (
    STORE / "gwas/ukbb_ppp/ppp_heritability/analysis/ppp_heritability_hapmap_3.parquet"
)
LOG_DIR = REPO / "experiments/claude/csf_ldsc/logs"

# Columns shared by both tables (the CSF table adds a precomputed p; we recompute it the
# same way for both so the comparison uses one code path).
H2_COL = "h2"
H2_SE_COL = "h2_se"
MEAN_CHI2_COL = "mean_chi2"
N_BAR_COL = "n_bar"
VARIANT_SET_COL = "variant_set"
ALL_VARIANTS = "all_variants"


@frozen
class Diagnostics:
    label: str
    n: int
    median_n_bar: float
    h2_mean: float
    h2_median: float
    h2_std: float
    h2_min: float
    h2_max: float
    frac_h2_negative: float
    frac_mean_chi2_below_1: float
    mean_chi2_mean: float
    corr_h2_mean_chi2: float
    n_p05: int
    n_p001: int
    n_bonferroni: int
    n_bonferroni_positive: int


def two_sided_wald_p(h2: np.ndarray, h2_se: np.ndarray) -> np.ndarray:
    """z = h2/se, p = 2*norm.sf(|z|): the same convention CsfProteinHeritabilityTask
    writes, recomputed here so PPP (which has no p column) is scored identically."""
    return 2.0 * scipy.stats.norm.sf(np.abs(h2 / h2_se))


def diagnostics(frame: pl.DataFrame, label: str) -> Diagnostics:
    h2 = frame[H2_COL].to_numpy()
    se = frame[H2_SE_COL].to_numpy()
    mc = frame[MEAN_CHI2_COL].to_numpy()
    n_bar = frame[N_BAR_COL].to_numpy()
    p = two_sided_wald_p(h2, se)
    bonferroni = 0.05 / len(frame)
    return Diagnostics(
        label=label,
        n=len(frame),
        median_n_bar=float(np.median(n_bar)),
        h2_mean=float(h2.mean()),
        h2_median=float(np.median(h2)),
        h2_std=float(h2.std()),
        h2_min=float(h2.min()),
        h2_max=float(h2.max()),
        frac_h2_negative=float((h2 < 0).mean()),
        frac_mean_chi2_below_1=float((mc < 1).mean()),
        mean_chi2_mean=float(mc.mean()),
        corr_h2_mean_chi2=float(np.corrcoef(h2, mc)[0, 1]),
        n_p05=int((p < 0.05).sum()),
        n_p001=int((p < 0.001).sum()),
        n_bonferroni=int((p < bonferroni).sum()),
        n_bonferroni_positive=int(((p < bonferroni) & (h2 > 0)).sum()),
    )


def all_variants_rows(frame: pl.DataFrame) -> pl.DataFrame:
    """The all-variants rows. CSF has only these; PPP also has cis_excluded rows (nearly
    identical numbers), which we drop so the comparison is like-for-like."""
    if VARIANT_SET_COL in frame.columns:
        return frame.filter(pl.col(VARIANT_SET_COL) == ALL_VARIANTS)
    return frame


def format_table(rows: list[Diagnostics]) -> str:
    labels = [d.label for d in rows]
    metrics: list[tuple[str, list[str]]] = [
        ("median N", [f"{d.median_n_bar:,.0f}" for d in rows]),
        ("n proteins/aptamers", [f"{d.n:,}" for d in rows]),
        ("h2 mean", [f"{d.h2_mean:.4f}" for d in rows]),
        ("h2 median", [f"{d.h2_median:.4f}" for d in rows]),
        ("h2 std", [f"{d.h2_std:.4f}" for d in rows]),
        ("h2 min", [f"{d.h2_min:.3f}" for d in rows]),
        ("h2 max", [f"{d.h2_max:.3f}" for d in rows]),
        ("fraction h2 < 0", [f"{d.frac_h2_negative:.3f}" for d in rows]),
        ("fraction mean_chi2 < 1", [f"{d.frac_mean_chi2_below_1:.3f}" for d in rows]),
        ("mean_chi2 (mean)", [f"{d.mean_chi2_mean:.4f}" for d in rows]),
        ("corr(h2, mean_chi2)", [f"{d.corr_h2_mean_chi2:.3f}" for d in rows]),
        ("n p<0.05", [f"{d.n_p05:,}" for d in rows]),
        ("  (chance)", [f"{0.05 * d.n:,.0f}" for d in rows]),
        ("n p<0.001", [f"{d.n_p001:,}" for d in rows]),
        ("  (chance)", [f"{0.001 * d.n:,.1f}" for d in rows]),
        ("survive Bonferroni", [f"{d.n_bonferroni:,}" for d in rows]),
        ("  of which h2 > 0", [f"{d.n_bonferroni_positive:,}" for d in rows]),
    ]
    name_w = max(len(m) for m, _ in metrics)
    col_w = max(12, *(len(x) for _, vals in metrics for x in vals))
    header = "statistic".ljust(name_w) + "".join(l.rjust(col_w) for l in labels)
    lines = [header, "-" * len(header)]
    for metric, vals in metrics:
        lines.append(metric.ljust(name_w) + "".join(v.rjust(col_w) for v in vals))
    return "\n".join(lines)


def main() -> None:
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    log = open(LOG_DIR / "h2_noise_floor_comparison.log", "w")

    def emit(*a: object) -> None:
        line = " ".join(str(x) for x in a)
        print(line)
        log.write(line + "\n")

    ppp = diagnostics(all_variants_rows(pl.read_parquet(PPP_PATH)), "UKB-PPP")
    csf = diagnostics(all_variants_rows(pl.read_parquet(CSF_PATH)), "CSF")

    emit("LDSC heritability diagnostics: UKB-PPP plasma vs Western CSF (all_variants)")
    emit("")
    emit(format_table([ppp, csf]))
    emit("")
    emit(
        "Reading: both tables' Bonferroni-significant hits are POSITIVE and PPP's few "
        "negatives hug zero -- negative h2 is the LDSC noise floor. CSF's ~10x smaller N "
        "puts mean_chi2 at the null, so the floor swallows ~half the estimates and none "
        "survive correction. It is a power problem, not a CSF-specific defect."
    )
    log.close()


if __name__ == "__main__":
    main()
