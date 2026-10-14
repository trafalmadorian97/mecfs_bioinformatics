"""How do MAF-bin annotations influence the polyfun prior?

Reads the ridge regression coefficients (gamma_raw) of the polyfun per-SNP
heritability prior on the baseline-LF functional annotations, isolates the 20
MAF-bin annotations (MAFbin_frequent_1..10 for common variants MAF>0.05, and
MAFbin_lowfreq_1..10 for low-frequency variants), and asks whether the
coefficients imply rarer variants are more likely causal.

Because each variant falls in exactly one MAF bin (verified here), the MAF-family
contribution to the linear predictor of a variant is simply that one bin's
gamma_raw, so comparing gamma_raw across bins compares the prior across the
frequency spectrum directly.

The frequent-vs-lowfreq split is definitional (lowfreq = rarer), so it answers
the rare/common question on its own. To order the bins WITHIN each group we join
one chromosome's annotation slice to the DecodeME genome-wide allele frequencies
and take the mean MAF per bin.

Run: pixi r python experiments/claude/maf_bin_prior/analyze_maf_bin_coefficients.py
"""

import polars as pl

from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

RIDGE = "assets/base_asset_store/reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annotation_ridge_weights/weights.parquet"
ANNOT_CHR1 = "assets/base_asset_store/reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annot_parquet_members/baselineLF2.2.UKB.1.annot.parquet"
# Build-37 (hg19), gwaslab-format harmonized dump, matching the annotation build.
DECODE = "assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_harmonized_dump_to_parquet.parquet"


def _bin_index(name: str) -> int:
    return int(name.rsplit("_", 1)[1])


def main() -> None:
    weights = pl.read_parquet(RIDGE)
    maf = weights.filter(pl.col("family") == "maf_bins").select(
        "annotation", "gamma_raw", "gamma_standardized"
    )
    freq_names = [n for n in maf["annotation"] if n.startswith("MAFbin_frequent_")]
    low_names = [n for n in maf["annotation"] if n.startswith("MAFbin_lowfreq_")]
    gamma = dict(zip(maf["annotation"], maf["gamma_raw"]))

    print("=== gamma_raw by MAF bin (common variants, MAF>0.05) ===")
    for n in sorted(freq_names, key=_bin_index):
        print(f"  {n:22s} gamma_raw={gamma[n]: .4e}")
    print("=== gamma_raw by MAF bin (low-frequency variants) ===")
    for n in sorted(low_names, key=_bin_index):
        print(f"  {n:22s} gamma_raw={gamma[n]: .4e}")

    fmean = sum(gamma[n] for n in freq_names) / len(freq_names)
    lmean = sum(gamma[n] for n in low_names) / len(low_names)
    print("\n=== group summary (definitional: lowfreq = rarer) ===")
    print(f"  mean gamma_raw, frequent (common): {fmean: .4e}")
    print(f"  mean gamma_raw, lowfreq  (rare)  : {lmean: .4e}")
    print(f"  lowfreq - frequent               : {lmean - fmean: .4e}")

    # Establish within-group MAF ordering empirically on chr1: which physical MAF
    # does each bin correspond to?
    maf_cols = freq_names + low_names
    annot = pl.read_parquet(ANNOT_CHR1, columns=["CHR", "BP", "A1", "A2", *maf_cols])
    per_variant_bins = annot.select(pl.sum_horizontal(maf_cols).alias("n_bins"))
    counts = per_variant_bins["n_bins"].value_counts().sort("n_bins")
    print("\n=== MAF-bin membership per variant (chr1 annotation slice) ===")
    print(counts)

    annot = annot.with_columns(
        unordered_allele_key("A1", "A2").alias("akey")
    ).rename({"BP": "POS"})
    decode = (
        pl.read_parquet(DECODE, columns=["CHR", "POS", "EA", "NEA", "EAF"])
        .filter(pl.col("CHR") == 1)
        .with_columns(
            unordered_allele_key("EA", "NEA").alias("akey"),
            pl.min_horizontal(pl.col("EAF"), 1.0 - pl.col("EAF")).alias("maf"),
        )
    )
    joined = annot.join(decode.select("CHR", "POS", "akey", "maf"), on=["CHR", "POS", "akey"], how="inner")
    print(f"\nchr1 annotation variants joined to DecodeME MAF: {joined.height}")
    rows = []
    for n in maf_cols:
        sub = joined.filter(pl.col(n) == 1)
        if sub.height:
            rows.append((n, sub.height, float(sub["maf"].mean()), gamma[n]))
    order = pl.DataFrame(
        rows, schema=["bin", "n", "mean_maf", "gamma_raw"], orient="row"
    ).sort("mean_maf")
    print("\n=== bins ordered by empirical mean MAF (rarest first) ===")
    with pl.Config(tbl_rows=40, fmt_float="full"):
        print(order)


if __name__ == "__main__":
    main()
