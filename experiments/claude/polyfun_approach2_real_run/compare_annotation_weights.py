"""Compare DecodeME Approach-2 annotation coefficients (tau_odd, tau_even) with the
Approach-1 ridge surrogate (ridge fit of the precomputed PolyFun snpvar_bin on the
same 187 baseline-LF annotations).

The raw scales differ (tau is per-SNP heritability; the Approach-1 gamma predicts
snpvar_bin), so coefficients are compared on the standardized scale
(gamma_standardized) by correlation, top-k overlap, and per-family sums. The
implied genome-wide priors are compared by joining Approach-1 snpvar_bin to
Approach-2 snpvar on the exact (CHR, POS, NEA, EA) key.

Run: pixi r python experiments/claude/polyfun_approach2_real_run/compare_annotation_weights.py \
    | tee experiments/claude/polyfun_approach2_real_run/compare_annotation_weights.log
"""

from pathlib import Path

import polars as pl
from scipy.stats import pearsonr, spearmanr

from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_STANDARDIZED_COL,
)
from mecfs_bio.constants.polyfun_constants import POLYFUN_TO_GWASLAB_KEY_RENAME

STORE = Path("assets/base_asset_store")
A1_WEIGHTS = (
    STORE
    / "reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annotation_ridge_weights/weights.parquet"
)
A2_DIR = (
    STORE
    / "gwas/ME_CFS/DecodeME/analysis/polyfun_l2_sldsc_prior/decode_me_gwas_1_l2_sldsc_snpvar"
)
A1_PRIOR = (
    STORE
    / "reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet"
)
KEY = ["CHR", "POS", "NEA", "EA"]
TOP_K = 20


def _std(path: Path, name: str) -> pl.DataFrame:
    return pl.read_parquet(path).select(
        ANNOTATION_COL, FAMILY_COL, pl.col(GAMMA_STANDARDIZED_COL).alias(name)
    )


def main() -> None:
    a1 = _std(A1_WEIGHTS, "approach1")
    odd = _std(A2_DIR / "weights_tau_odd.parquet", "tau_odd").drop(FAMILY_COL)
    even = _std(A2_DIR / "weights_tau_even.parquet", "tau_even").drop(FAMILY_COL)
    w = a1.join(odd, on=ANNOTATION_COL).join(even, on=ANNOTATION_COL)
    assert w.height == a1.height == odd.height == even.height, "annotation sets differ"
    print(f"annotations compared: {w.height}\n")

    print("== standardized-coefficient correlation (Pearson / Spearman) ==")
    pairs = [("tau_odd", "tau_even"), ("approach1", "tau_odd"), ("approach1", "tau_even")]
    for a, b in pairs:
        p = pearsonr(w[a], w[b]).statistic
        s = spearmanr(w[a], w[b]).statistic
        print(f"  {a:>9} vs {b:<9}  pearson={p:+.3f}  spearman={s:+.3f}")

    print(f"\n== top-{TOP_K} annotations by standardized coefficient: overlap ==")
    tops = {c: set(w.sort(c, descending=True).head(TOP_K)[ANNOTATION_COL]) for c in ("approach1", "tau_odd", "tau_even")}
    for a, b in pairs:
        print(f"  {a:>9} & {b:<9}  {len(tops[a] & tops[b])}/{TOP_K}")

    print("\n== sign agreement ==")
    for a, b in pairs:
        agree = (w[a].sign() == w[b].sign()).mean()
        print(f"  {a:>9} vs {b:<9}  {agree:.2f}")

    print("\n== per-family sum of standardized coefficients (each column scaled by its max |value|) ==")
    fam = w.group_by(FAMILY_COL).agg(pl.col("approach1", "tau_odd", "tau_even").sum())
    fam = fam.with_columns(
        (pl.col(c) / pl.col(c).abs().max()).round(2) for c in ("approach1", "tau_odd", "tau_even")
    ).sort("approach1", descending=True)
    print(fam)

    for c in ("approach1", "tau_odd", "tau_even"):
        print(f"\n== top 8 / bottom 4 annotations: {c} ==")
        s = w.sort(c, descending=True).select(ANNOTATION_COL, FAMILY_COL, c)
        print(s.head(8))
        print(s.tail(4))

    print("\n== implied genome-wide prior: Approach-1 snpvar_bin vs Approach-2 snpvar ==")
    a1_prior = pl.scan_parquet(A1_PRIOR).rename(POLYFUN_TO_GWASLAB_KEY_RENAME).select(*KEY, "snpvar_bin")
    a2_prior = pl.scan_parquet(A2_DIR / "snpvar.parquet").select(*KEY, "snpvar")
    joined = a1_prior.join(a2_prior, on=KEY, how="inner").collect()
    print(f"  variants in both: {joined.height:,}")
    print(
        f"  pearson={pearsonr(joined['snpvar_bin'], joined['snpvar']).statistic:+.3f}"
        f"  spearman={spearmanr(joined['snpvar_bin'], joined['snpvar']).statistic:+.3f}"
    )
    for chrom_parity, label in ((1, "odd chrs (scored by tau_even)"), (0, "even chrs (scored by tau_odd)")):
        sub = joined.filter(pl.col("CHR") % 2 == chrom_parity)
        print(f"  {label}: spearman={spearmanr(sub['snpvar_bin'], sub['snpvar']).statistic:+.3f}")


if __name__ == "__main__":
    main()
