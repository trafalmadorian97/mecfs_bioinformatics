"""Net effect on the SUSIE inputs: for each current DecodeME polyfun-explain harmonized locus file,
compare every row's final EA/NEA with the liftover-37 table (pre gwaslab/ANNOVAR/harmonizer) by
SNPID. Rows whose indel alleles end up swapped relative to the source are insertion<->deletion
relabels that the final harmonizer did NOT undo. Also report whether those rows are in any
credible set of the corresponding SUSIE runs."""
from pathlib import Path

import polars as pl

D = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME")
src = (
    pl.scan_parquet(D / "processed/decode_me_gwas_1_liftover_to_37_parquet_file.parquet")
    .select("SNPID", pl.col("EA").cast(pl.Utf8).alias("s_EA"), pl.col("NEA").cast(pl.Utf8).alias("s_NEA"))
    .collect()
)
for p in sorted((D / "processed").glob("decode_me_polyfun_explain*harmonized_with_ref.parquet")):
    h = pl.read_parquet(p).join(src, on="SNPID")
    indel = (h["EA"].str.len_chars() > 1) | (h["NEA"].str.len_chars() > 1)
    swapped = (h["EA"] == h["s_NEA"]) & (h["NEA"] == h["s_EA"])
    stem = p.name.removesuffix("_gwas_harmonized_with_ref.parquet")
    relabeled = h.filter(swapped & indel).select("SNPID", "POS", "s_EA", "s_NEA", "EA", "NEA")
    in_cs = []
    for run in sorted((D / "analysis").glob(stem + "_*_susie_*")):
        cs = run / "combined_cs.parquet"
        if cs.exists() and relabeled.height:
            hit = pl.read_parquet(cs).with_columns(pl.col("POS").cast(pl.Int64), pl.col("EA").cast(pl.Utf8), pl.col("NEA").cast(pl.Utf8)).join(relabeled.select("POS", "EA", "NEA"), on=["POS", "EA", "NEA"], how="semi")
            if hit.height:
                in_cs.append((run.name, hit.height))
    print(f"{stem}: rows={h.height:,} indels={int(indel.sum()):,} net-swapped indels={relabeled.height} net-swapped SNVs={int((swapped & ~indel).sum())} in credible sets={in_cs}")
    if relabeled.height:
        print(relabeled)
