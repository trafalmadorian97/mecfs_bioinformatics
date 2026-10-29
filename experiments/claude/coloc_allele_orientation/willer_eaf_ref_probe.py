"""
Probe: what is the eaf_ref column of the Willer et al. 2013 LDL GWAS (GCST002222 harmonised)?

Hypothesis: eaf_ref is GLGC's Freq.A1.1000G.EUR, the 1000 Genomes EUR frequency of the
ORIGINAL effect_allele (A1). Test by joining non-palindromic biallelic SNVs to the repo's
1000G EUR hg38 allele-frequency panel at the GWAS Catalog GRCh38 position (hm_chrom, hm_pos) and comparing eaf_ref
with the panel frequency of effect_allele vs the panel frequency of other_allele. Palindromes
are excluded because their allele identity cannot be read off the panel's REF/ALT.

Also reports whether beta is relative to effect_allele, via the GWAS Catalog's own sign
relationship between beta and hm_beta when its hm_effect_allele equals effect_allele.

Usage:
  pixi r python experiments/claude/coloc_allele_orientation/willer_eaf_ref_probe.py \
      2>&1 | tee experiments/claude/logs/willer_eaf_ref_probe.log
"""

import polars as pl

WILLER = (
    "assets/base_asset_store/gwas/ldl/willer_et_al_2013/raw/"
    "24097068-GCST002222-EFO_0004611.h.tsv.gz"
)
PANEL = (
    "assets/base_asset_store/reference_data/reference_panel_allele_frequencies/"
    "eur_hg38_30x/processed/panel_allele_frequencies.parquet"
)
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def main() -> None:
    willer = (
        pl.scan_csv(WILLER, separator="\t", null_values=["NA"], infer_schema=False)
        .select(
            chrom=pl.col("hm_chrom").cast(pl.Int32, strict=False),
            pos=pl.col("hm_pos").cast(pl.Int32, strict=False),
            ea=pl.col("effect_allele").str.to_uppercase(),
            oa=pl.col("other_allele").str.to_uppercase(),
            eaf_ref=pl.col("eaf_ref").cast(pl.Float64),
            beta=pl.col("beta").cast(pl.Float64),
            hm_beta=pl.col("hm_beta").cast(pl.Float64),
            hm_ea=pl.col("hm_effect_allele"),
        )
        .filter(pl.col("ea").str.len_bytes() == 1, pl.col("oa").str.len_bytes() == 1)
        .filter(pl.col("ea").replace_strict(COMPLEMENT, default=None) != pl.col("oa"))
        .collect()
    )
    print(f"non-palindromic SNVs: {willer.height:,}")
    print(f"eaf_ref range: {willer['eaf_ref'].min()} .. {willer['eaf_ref'].max()}")
    print(f"fraction eaf_ref <= 0.5: {(willer['eaf_ref'] <= 0.5).mean():.4f}")

    panel = pl.scan_parquet(PANEL).filter(
        pl.col("REF").str.len_bytes() == 1, pl.col("ALT").str.len_bytes() == 1
    )
    joined = (
        willer.lazy()
        .join(panel, left_on=["chrom", "pos"], right_on=["CHR", "POS"], how="inner")
        .with_columns(
            # Panel frequency of each Willer allele; strand-flipped rows allowed.
            ea_freq=pl.when(pl.col("ea") == pl.col("ALT"))
            .then(pl.col("AF"))
            .when(pl.col("ea") == pl.col("REF"))
            .then(1 - pl.col("AF"))
            .when(pl.col("ea").replace_strict(COMPLEMENT, default=None) == pl.col("ALT"))
            .then(pl.col("AF"))
            .when(pl.col("ea").replace_strict(COMPLEMENT, default=None) == pl.col("REF"))
            .then(1 - pl.col("AF")),
        )
        .filter(pl.col("ea_freq").is_not_null())
        .unique(["chrom", "pos", "ea", "oa"], keep="none")
        .collect()
    )
    print(f"matched to 1000G EUR hg38 panel : {joined.height:,}")
    d_ea = (joined["eaf_ref"] - joined["ea_freq"]).abs()
    d_oa = (joined["eaf_ref"] - (1 - joined["ea_freq"])).abs()
    print("eaf_ref vs panel freq of effect_allele:")
    print(
        f"  corr={pl.DataFrame({'a': joined['eaf_ref'], 'b': joined['ea_freq']}).select(pl.corr('a', 'b')).item():.4f}"
        f"  median|diff|={d_ea.median():.4f}  frac|diff|<0.05={(d_ea < 0.05).mean():.4f}"
    )
    print("eaf_ref vs panel freq of other_allele:")
    print(
        f"  median|diff|={d_oa.median():.4f}  frac|diff|<0.05={(d_oa < 0.05).mean():.4f}"
    )
    print(f"fraction closer to effect_allele: {(d_ea < d_oa).mean():.4f}")

    same = willer.filter(pl.col("hm_ea") == pl.col("ea"), pl.col("hm_beta").is_not_null())
    flipped = willer.filter(
        pl.col("hm_ea") == pl.col("oa"), pl.col("hm_beta").is_not_null()
    )
    print(
        "GWAS Catalog sign check: hm_beta == beta when hm_ea == ea: "
        f"{(same['hm_beta'] == same['beta']).mean():.4f} (n={same.height:,}); "
        f"hm_beta == -beta when hm_ea == oa: "
        f"{(flipped['hm_beta'] == -flipped['beta']).mean():.4f} (n={flipped.height:,})"
    )


if __name__ == "__main__":
    main()
