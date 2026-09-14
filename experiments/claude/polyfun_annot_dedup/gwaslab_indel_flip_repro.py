"""Minimal reproduction: gwaslab harmonization flips a deletion into an insertion at a repeat.

A GWAS reports a deletion at hg19 chr15:54698192 (reference TGGGGGG): REF=TG, ALT=T, with the
deletion allele T as the effect allele (EA=T, NEA=TG, EAF=0.845). The 1000 Genomes EUR VCF used for
strand inference contains BOTH
    REF=T,  ALT=TG, AF=0.0    (an insertion of G; monomorphic in EUR)
    REF=TG, ALT=T,  AF=0.862  (the deletion; matches the GWAS frequency)
Because both alleles T and TG match the reference, check_ref marks the variant "indistinguishable"
(STATUS digit 6 = 6). infer_strand's check_unkonwn_indel then returns "flip" (digit 7 = 6) on the
FIRST VCF record whose alleles are the string-swap of the sumstats alleles within a 0.2 frequency
tolerance -- the AF=0 insertion -- and _flip_allele_stats swaps EA/NEA and negates BETA.
The matching deletion record is never considered.

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/gwaslab_indel_flip_repro.py
"""

from pathlib import Path

import gwaslab as gl
import pandas as pd
import pysam

FASTA = Path.home() / ".gwaslab/hg19.fa"
VCF = Path.home() / ".gwaslab/EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"  # gl.get_path("1kg_eur_hg19")
CHROM, POS = 15, 54698192
EA, NEA, EAF, BETA, SE = "T", "TG", 0.845, -0.035, 0.0175

print(f"gwaslab version: {gl.__version__}")

# 1. Reference sequence at the site.
ref_seq = pysam.FastaFile(str(FASTA)).fetch(f"chr{CHROM}", POS - 1, POS + 6).upper()
print(f"\n1. hg19 chr{CHROM}:{POS}..{POS + 6} = {ref_seq}")
assert ref_seq.startswith(EA) and ref_seq.startswith(NEA), "both alleles should match the reference"

# 2. VCF records at the site.
records = {
    (rec.ref, alt): rec.info["AF"][i]
    for rec in pysam.VariantFile(str(VCF)).fetch(str(CHROM), POS - 1, POS)
    if rec.pos == POS
    for i, alt in enumerate(rec.alts)
}
print(f"\n2. 1kg EUR VCF records at {CHROM}:{POS} (REF, ALT) -> AF: {records}")
assert records[("T", "TG")] == 0.0, "expected monomorphic insertion REF=T ALT=TG AF=0.0"
assert abs(records[("TG", "T")] - 0.862) < 0.001, "expected deletion REF=TG ALT=T AF~0.862"

# 3. One-row sumstats through gwaslab's harmonization pipeline.
sumstats = gl.Sumstats(
    pd.DataFrame(
        {
            "SNPID": [f"{CHROM}:{POS}:{NEA}:{EA}"],
            "CHR": [CHROM],
            "POS": [POS],
            "EA": [EA],
            "NEA": [NEA],
            "EAF": [EAF],
            "BETA": [BETA],
            "SE": [SE],
            "N": [275_488],
        }
    ),
    snpid="SNPID",
    chrom="CHR",
    pos="POS",
    ea="EA",
    nea="NEA",
    eaf="EAF",
    beta="BETA",
    se="SE",
    n="N",
    build="19",
    verbose=False,
)
sumstats.harmonize(
    basic_check=True,
    ref_seq=str(FASTA),
    ref_infer=str(VCF),
    ref_alt_freq="AF",
    threads=1,
    verbose=False,
)
out = sumstats.data.iloc[0]
print("\n3. Before gwaslab harmonization: EA=T  NEA=TG  EAF=0.845  BETA=-0.035")
print(
    f"   After  gwaslab harmonization: EA={out['EA']}  NEA={out['NEA']}  "
    f"EAF={out['EAF']:.3f}  BETA={out['BETA']:.3f}  STATUS={out['STATUS']}"
)
flipped = out["EA"] == NEA and out["NEA"] == EA
print(
    "\n   RESULT: gwaslab flipped the variant to REF=T ALT=TG (the AF=0.0 insertion record)."
    if flipped
    else "\n   RESULT: not flipped (bug not reproduced)."
)
assert flipped
