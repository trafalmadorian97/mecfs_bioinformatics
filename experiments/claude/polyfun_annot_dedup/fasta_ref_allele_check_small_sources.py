"""Convention check (which allele column matches hg19) over ALL rows of the small sources that have
no chr22 data: the local UKBB LD label windows and the DecodeME harmonized polyfun-explain loci."""
from collections import Counter

import polars as pl

from experiments.claude.polyfun_annot_dedup.fasta_ref_allele_check import (
    FASTA,
    Genome,
    _matches,
    decode_harmonized,
    ld_labels,
)

genome = Genome(FASTA["hg19"])
for name, frame, col_a, col_b in [
    ("UKBB Broad LD labels (hg19)", ld_labels(), "allele1", "allele2"),
    ("DecodeME harmonized polyfun-explain loci (hg19)", decode_harmonized(), "NEA", "EA"),
]:
    print(f"\n=== {name}: a={col_a}, b={col_b} ===")
    counts: dict[str, Counter[str]] = {"SNV": Counter(), "indel": Counter()}
    for chrom, pos, a, b in frame.sort("chr", "pos").select("chr", "pos", "a", "b").iter_rows():
        seq = genome.seq(chrom)
        kind = "SNV" if len(a) == 1 and len(b) == 1 else "indel"
        ma, mb = _matches(seq, pos, a), _matches(seq, pos, b)
        counts[kind]["both" if ma and mb else f"{col_a}_only" if ma else f"{col_b}_only" if mb else "neither"] += 1
    for kind, c in counts.items():
        total = sum(c.values())
        print(f"  {kind} (n={total:,}): " + ", ".join(f"{k}={v:,} ({100 * v / total:.2f}%)" for k, v in c.most_common()))
