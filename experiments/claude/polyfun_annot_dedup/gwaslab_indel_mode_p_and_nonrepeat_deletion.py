"""Two questions bearing on the CleanIndelOrientation design.

A. Does infer_strand_kwargs={"mode": "p"} stop gwaslab harmonization from flipping the chr15:54698192
   deletion (see gwaslab_indel_flip_repro.py)?
B. Does check_ref mark an ordinary NON-repeat deletion (REF=XY, ALT=X, next reference base != Y) as
   "both alleles on genome, indistinguishable" (STATUS digit 6 = 6)? If so, "both alleles match the
   FASTA" is true of every deletion, not only of repeat sites.

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/gwaslab_indel_mode_p_and_nonrepeat_deletion.py \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/gwaslab_indel_mode_p_and_nonrepeat_deletion.log
"""

from pathlib import Path

import gwaslab as gl
import pandas as pd
import pysam

FASTA = Path.home() / ".gwaslab/hg19.fa"
VCF = Path.home() / ".gwaslab/EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
CHROM = 15
REPEAT_POS = 54698192

fasta = pysam.FastaFile(str(FASTA))
window_start = REPEAT_POS + 1000
window = fasta.fetch(f"chr{CHROM}", window_start - 1, window_start + 200).upper()
offset = next(
    i for i in range(len(window) - 3) if len({window[i], window[i + 1], window[i + 2]}) == 3 and "N" not in window[i : i + 3]
)
del_pos = window_start + offset
del_ref, del_alt = window[offset : offset + 2], window[offset]
print(f"Non-repeat deletion site chr{CHROM}:{del_pos} context={window[offset : offset + 6]} REF={del_ref} ALT={del_alt}")


def harmonize(rows: list[dict], infer_strand_kwargs: dict) -> pd.DataFrame:
    df = pd.DataFrame(rows)
    df["SNPID"] = [f"{CHROM}:{r['POS']}:{r['NEA']}:{r['EA']}" for r in rows]
    df["CHR"] = CHROM
    df["SE"] = 0.0175
    df["N"] = 275_488
    s = gl.Sumstats(
        df, snpid="SNPID", chrom="CHR", pos="POS", ea="EA", nea="NEA", eaf="EAF", beta="BETA", se="SE", n="N",
        build="19", verbose=False,
    )
    s.harmonize(
        basic_check=True, ref_seq=str(FASTA), ref_infer=str(VCF), ref_alt_freq="AF", threads=1, verbose=False,
        infer_strand_kwargs=infer_strand_kwargs,
    )
    return s.data[["POS", "EA", "NEA", "EAF", "BETA", "STATUS"]]


rows = [
    {"POS": REPEAT_POS, "EA": "T", "NEA": "TG", "EAF": 0.845, "BETA": -0.035},
    {"POS": del_pos, "EA": del_alt, "NEA": del_ref, "EAF": 0.3, "BETA": 0.01},
]
print("\nInput:")
print(pd.DataFrame(rows))
for kwargs in [{}, {"mode": "p"}]:
    print(f"\ninfer_strand_kwargs={kwargs}:")
    out = harmonize(rows, kwargs)
    out["digit6"] = out["STATUS"].astype(str).str[5]
    out["digit7"] = out["STATUS"].astype(str).str[6]
    print(out.to_string())
