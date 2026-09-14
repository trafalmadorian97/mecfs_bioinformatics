"""Trace the single variant flipped by HarmonizeGWASWithReferenceViaAlleles at the DecodeME chr15
locus (SNPID 15:54405994:TG:T, hg19 15:54698192) through every DecodeME processing stage, and
compare with the UKBB LD panel, PolyFun prior/annotations, and UKB-PPP rows at the same site."""
from pathlib import Path

import polars as pl
import pysam

S = Path("assets/base_asset_store")
D = S / "gwas/ME_CFS/DecodeME"
HG38_POS, HG19_POS = 54405994, 54698192
cfg = pl.Config(tbl_cols=-1, tbl_width_chars=250)
cfg.__enter__()

print("hg38 ref at 15:54405994..+6:", pysam.FastaFile(str(Path.home() / ".gwaslab/hg38.fa")).fetch("chr15", HG38_POS - 1, HG38_POS + 6).upper())
print("hg19 ref at 15:54698192..+6:", pysam.FastaFile(str(Path.home() / ".gwaslab/hg19.fa")).fetch("chr15", HG19_POS - 1, HG19_POS + 6).upper())

print("\n-- DecodeME raw regenie (hg38) --")
print(
    pl.scan_csv(D / "raw/DecodeME Summary Statistics/gwas_1.regenie.gz", separator=" ", schema_overrides={"CHROM": pl.Utf8})
    .filter((pl.col("CHROM") == "15") & pl.col("GENPOS").is_between(HG38_POS - 2, HG38_POS + 2))
    .select("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "BETA", "SE")
    .collect()
)
for name in [
    "DecodeME_keep_version_parquet_table_from_sumstats.parquet",
    "decode_me_gwas_1_liftover_to_37_parquet_file.parquet",
    "decode_me_gwas_1_assign_rsids_via_dbsnp150.parquet",
    "decode_me_polyfun_explainchr15_54500000_55500000_palindromes_keep_gwas_harmonized_with_ref.parquet",
]:
    lf = pl.scan_parquet(D / "processed" / name)
    cols = [c for c in ["SNPID", "rsID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "REF", "ALT"] if c in lf.collect_schema().names()]
    print(f"\n-- {name} --")
    print(lf.filter(pl.col("SNPID").str.starts_with("15:54405994:")).select(cols).with_columns(pl.col(pl.Categorical).cast(pl.Utf8)).collect())

print("\n-- UKBB Broad LD labels (hg19) --")
print(pl.read_csv(S / "reference_data/ukbb_reference_ld/chr15_53000001_56000001/raw/chr15_53000001_56000001.gz", separator="\t").filter(pl.col("position").is_between(HG19_POS - 2, HG19_POS + 2)))

print("\n-- PolyFun prior (hg19) --")
print(pl.scan_parquet(S / "reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet").filter((pl.col("CHR") == 15) & pl.col("BP").is_between(HG19_POS - 2, HG19_POS + 2)).collect())

annot = pl.scan_parquet(S / "reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annot_parquet_members/baselineLF2.2.UKB.15.annot.parquet")
maf = [c for c in annot.collect_schema().names() if c.startswith("MAFbin")]
print("\n-- PolyFun annotations (hg19): MAF bins set to 1 --")
a = annot.filter(pl.col("BP").is_between(HG19_POS - 2, HG19_POS + 2)).select("BP", "SNP", "A1", "A2", *maf).collect()
for row in a.iter_rows(named=True):
    print(row["BP"], row["SNP"], row["A1"], row["A2"], [c for c in maf if row[c] == 1])

print("\n-- UKB-PPP RABGAP1L (ID carries hg19 pos) --")
print(
    pl.scan_parquet(S / "reference_data/ukbb_ppp/sun_et_al_2023/extracted/ukbb_ppp_rabgap1l_sumstats_stacked.parquet")
    .filter(pl.col("CHROM") == 15)
    .filter(pl.col("ID").str.split(":").list.get(1).cast(pl.Int64).is_between(HG19_POS - 2, HG19_POS + 2))
    .select("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "INFO")
    .collect()
)
