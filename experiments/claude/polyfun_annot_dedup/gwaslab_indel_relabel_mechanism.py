"""Why does gwaslab harmonization swap EA/NEA of some DecodeME indels?

gwaslab 4.2.2 mechanism (hm/hm_harmonize_sumstats.py):
  1. check_ref: an indel whose EA and NEA both match the FASTA at POS gets STATUS digit 6 = 6
     ("both alleles on genome but indistinguishable").
  2. infer_strand (mode "pi", "i" part): for those indels, check_unkonwn_indel scans 1kg EUR VCF
     records overlapping POS and returns digit 7 =
        3 (keep)  if a record has pos==POS, REF==NEA, EA in ALTs and |AF - EAF| < 0.2
        6 (flip)  if a record has pos==POS, REF==EA, NEA in ALTs and |AF - (1-EAF)| < 0.2
        8         if any overlapping record has MAF > 0.4, or nothing matched.
  3. flip_allele_stats: digit 5 in [1,2,3], digit 6 in [6,7], digit 7 == 6 -> swap EA/NEA,
     negate BETA, EAF -> 1-EAF; digit 7 becomes 4.

The flip branch equates allele STRINGS across two different variant records: sumstats "EA=T,
NEA=TG" (deletion of G) is treated as the same variant as VCF "REF=T, ALT=TG" (insertion of G).
At a repeat these are different haplotypes. The only guard is a 0.2 allele-frequency tolerance.

This script, for every indel gwaslab harmonization swapped (liftover-37 vs harmonized dump):
  - re-runs gwaslab's own check_unkonwn_indel with the current 1kg VCF;
  - lists the VCF records at POS: the "opposite-representation" record gwaslab matched, and whether
    the VCF ALSO has a same-representation record (the one that describes DecodeME's variant);
  - compares DecodeME's original EAF with each record's AF.
It also tabulates STATUS digit 7 over all indistinguishable indels.

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/gwaslab_indel_relabel_mechanism.py \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/gwaslab_indel_relabel_mechanism.log
"""

from pathlib import Path

import polars as pl
from gwaslab.hm.hm_harmonize_sumstats import check_unkonwn_indel
from pysam import VariantFile

D = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME/processed")
VCF = Path.home() / ".gwaslab/EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
AF = "AF"


def load(name: str) -> pl.DataFrame:
    return (
        pl.scan_parquet(D / name)
        .select("SNPID", "CHR", "POS", "EA", "NEA", "EAF", "BETA", "STATUS")
        .with_columns(pl.col(pl.Categorical).cast(pl.Utf8), pl.col("STATUS").cast(pl.Utf8))
        .collect()
    )


pre = load("decode_me_gwas_1_liftover_to_37_parquet_file.parquet")
post = load("decode_me_gwas_1_harmonized_dump_to_parquet.parquet")
j = post.join(
    pre.select("SNPID", *[pl.col(c).alias(f"pre_{c}") for c in ["EA", "NEA", "EAF", "BETA", "STATUS"]]),
    on="SNPID",
)
indel = (pl.col("pre_EA").str.len_chars() > 1) | (pl.col("pre_NEA").str.len_chars() > 1)
j = j.filter(indel).with_columns(
    pl.col("STATUS").str.slice(5, 1).alias("d6"),
    pl.col("STATUS").str.slice(6, 1).alias("d7"),
    ((pl.col("EA") == pl.col("pre_NEA")) & (pl.col("NEA") == pl.col("pre_EA"))).alias("swapped"),
)
print("=== All indels: STATUS digit 6 (check_ref) x digit 7 (infer_strand) x swapped ===")
with pl.Config(tbl_rows=40):
    print(j.group_by("d6", "d7", "swapped").len().sort("d6", "d7", "swapped"))

sw = j.filter(pl.col("swapped"))
print(f"\n=== Swapped indels: {sw.height:,} ===")
vcf = VariantFile(str(VCF))
rows = []
for r in sw.iter_rows(named=True):
    chrom, pos = str(r["CHR"]), int(r["POS"])
    ea, nea, eaf = r["pre_EA"], r["pre_NEA"], float(r["pre_EAF"])
    # gwaslab passes ref=NEA, alt=EA, start=POS-1, end=POS, status as before infer_strand (digit 7 = 9).
    status_in = int(r["pre_STATUS"][:5] + "69")
    new_status = check_unkonwn_indel(chrom, pos - 1, pos, nea, ea, eaf, vcf, AF, 0.4, status_in, None, 0.2)
    opposite_af = same_af = None
    n_records = 0
    for rec in vcf.fetch(chrom, pos - 1, pos):
        n_records += 1
        if rec.pos != pos:
            continue
        if rec.ref == ea and nea in rec.alts:
            opposite_af = rec.info[AF][rec.alts.index(nea)]
        if rec.ref == nea and ea in rec.alts:
            same_af = rec.info[AF][rec.alts.index(ea)]
    rows.append(
        {
            "SNPID": r["SNPID"],
            "d6": r["d6"],
            "d7": r["d7"],
            "POS": pos,
            "pre_EA": ea,
            "pre_NEA": nea,
            "pre_EAF": eaf,
            "gwaslab_digit7_now": str(new_status)[6],
            "vcf_opposite_record_AF_of_NEA": opposite_af,
            "vcf_same_record_AF_of_EA": same_af,
            "n_overlapping_records": n_records,
        }
    )
res = pl.DataFrame(rows).with_columns(
    (pl.col("vcf_opposite_record_AF_of_NEA") - (1 - pl.col("pre_EAF"))).abs().alias("daf_opposite"),
    (pl.col("vcf_same_record_AF_of_EA") - pl.col("pre_EAF")).abs().alias("daf_same"),
    ((pl.col("pre_EA").str.len_chars() > pl.col("pre_NEA").str.len_chars())).alias("pre_is_insertion"),
)
print("gwaslab check_unkonwn_indel digit 7 with CURRENT VCF (6 = flip reproduced):")
print(res.group_by("gwaslab_digit7_now").len().sort("gwaslab_digit7_now"))
print("\nVCF has opposite-representation record at POS:", int(res["vcf_opposite_record_AF_of_NEA"].is_not_null().sum()))
print("VCF ALSO has same-representation record at POS:", int(res["vcf_same_record_AF_of_EA"].is_not_null().sum()))
both = res.filter(pl.col("vcf_same_record_AF_of_EA").is_not_null() & pl.col("vcf_opposite_record_AF_of_NEA").is_not_null())
print(
    f"  of those with both records, the same-representation record fits DecodeME's frequency better: "
    f"{int((both['daf_same'] < both['daf_opposite']).sum())}/{both.height}"
)
print("\nDAF of the opposite record that triggered the flip (tolerance 0.2):")
print(res["daf_opposite"].describe())
print("\nSwapped indels by original type (pre_is_insertion):")
print(res.group_by("pre_is_insertion").len())
with pl.Config(tbl_cols=-1, tbl_width_chars=250, tbl_rows=20):
    print(res.filter(pl.col("SNPID").is_in(["15:54405994:TG:T", "1:7857852:G:GA", "15:55098921:CA:C", "17:52731691:C:CA"])))
    print(res.sort("daf_opposite", descending=True, nulls_last=True).head(8))

print("\n=== Restricted to the indistinguishable-indel branch (digit6 == 6, digit7 == 4) ===")
ib = res.filter((pl.col("d6") == "6") & (pl.col("d7") == "4"))
first_match_opposite = []
for r in ib.iter_rows(named=True):
    # Which record does gwaslab's loop hit first, and is it the opposite representation?
    for rec in vcf.fetch(str(sw.filter(pl.col("SNPID") == r["SNPID"])["CHR"][0]), r["POS"] - 1, r["POS"]):
        if rec.pos == r["POS"] and ((rec.ref == r["pre_NEA"] and r["pre_EA"] in rec.alts) or (rec.ref == r["pre_EA"] and r["pre_NEA"] in rec.alts)):
            first_match_opposite.append(rec.ref == r["pre_EA"])
            break
    else:
        first_match_opposite.append(None)
ib = ib.with_columns(pl.Series("first_matching_record_is_opposite", first_match_opposite))
print(f"n = {ib.height:,}")
print("flip reproduced with current VCF:", ib.group_by("gwaslab_digit7_now").len().sort("gwaslab_digit7_now").rows())
print("pre-gwaslab type (insertion=True):", ib.group_by("pre_is_insertion").len().rows())
print("VCF also has the same-representation record:", int(ib["vcf_same_record_AF_of_EA"].is_not_null().sum()))
b2 = ib.filter(pl.col("vcf_same_record_AF_of_EA").is_not_null() & pl.col("vcf_opposite_record_AF_of_NEA").is_not_null())
print(f"  same-representation record fits DecodeME EAF better: {int((b2['daf_same'] < b2['daf_opposite']).sum())}/{b2.height}")
print("first matching record in gwaslab's scan is the opposite one:", ib.group_by("first_matching_record_is_opposite").len().rows())
print("opposite record AF == 0 (monomorphic in 1kg EUR):", int((ib["vcf_opposite_record_AF_of_NEA"] == 0).sum()))
print("DAF of opposite record:", ib["daf_opposite"].describe().rows())
