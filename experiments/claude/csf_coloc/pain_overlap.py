"""Overlap between the fibromyalgia and multisite pain CSF coloc runs: top pairs per trait,
proteins tested in both, and overlapping trait loci.
Usage: pixi r python experiments/claude/csf_coloc/pain_overlap.py"""

import polars as pl

from mecfs_bio.analysis import csf_coloc_analysis as a
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.constants.protein_coloc_constants import (
    CIS_STATUS_COL,
    GENE_COL,
    LOCUS_ID_COL,
    PP_H3_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    STATUS_COL,
)

pl.Config.set_tbl_rows(60)
pl.Config.set_fmt_str_lengths(40)
COLS = [LOCUS_ID_COL, PROTEIN_ID_COL, GENE_COL, CIS_STATUS_COL, PP_H3_COL, PP_H4_COL]


def load(group) -> pl.DataFrame:
    out = DEFAULT_RUNNER.meta_to_path(group.coloc_task.meta)
    return pl.read_parquet(out / RAW_RESULTS_FILENAME)


def loci(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.select(LOCUS_ID_COL)
        .unique()
        .with_columns(
            pl.col(LOCUS_ID_COL)
            .str.extract_groups(r"^chr(\w+):(\d+)-(\d+)$")
            .struct.rename_fields(["chrom", "start", "end"])
            .alias("g")
        )
        .unnest("g")
        .with_columns(pl.col("start").cast(pl.Int64), pl.col("end").cast(pl.Int64))
    )


fibro = load(a.KERREBIJN_FIBRO_CSF_COLOC)
msp = load(a.JOHNSTON_PAIN_CSF_COLOC)
for name, df in (("fibro", fibro), ("msp", msp)):
    print(f"===== {name}: top tested pairs by PP.H4")
    print(df.filter(pl.col(STATUS_COL) == "tested").sort(PP_H4_COL, descending=True).select(COLS).head(10))

f_t = fibro.filter(pl.col(STATUS_COL) == "tested")
m_t = msp.filter(pl.col(STATUS_COL) == "tested")
shared = f_t.join(m_t, on=PROTEIN_ID_COL, suffix="_msp")
print("===== proteins tested for both traits")
print(shared.select(PROTEIN_ID_COL, GENE_COL, LOCUS_ID_COL, PP_H4_COL, LOCUS_ID_COL + "_msp", PP_H4_COL + "_msp"))

lf, lm = loci(fibro), loci(msp)
ov = lf.join(lm, on="chrom", suffix="_msp").filter(
    (pl.col("start") <= pl.col("end_msp")) & (pl.col("start_msp") <= pl.col("end"))
)
print(f"===== overlapping trait loci (fibro {lf.height}, msp {lm.height} loci)")
print(ov.select(LOCUS_ID_COL, LOCUS_ID_COL + "_msp"))
