"""
Which DecodeME variants did the switch to the gnomAD v2.1.1 nfe_nwe panel add, drop, or
re-orient in the harmonized tables?

Compares each baseline harmonized table (1000 Genomes EUR panel) with the rebuilt one, keyed
on SNPID, and classifies the variants as palindromic SNV, other SNV, or indel:

- rows only in the old table, and rows only in the new table;
- shared rows whose BETA sign changed. A palindrome whose strand decision reversed keeps its
  alleles after REF orientation but gets -BETA and 1-EAF; any other sign change comes with an
  effect-allele change.

Run after rebuild_decode_me.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.harmonized_row_changes \
        2>&1 | tee experiments/claude/decode_me_panel_switch/harmonized_row_changes.log
"""

import polars as pl

from experiments.claude.decode_me_panel_switch.compare_decode_me_to_baseline import (
    Record,
    load_records,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_SNPID_COL,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
COLUMNS = [GWASLAB_SNPID_COL, EA, NEA, GWASLAB_BETA_COL, GWASLAB_EFFECT_ALLELE_FREQ_COL]


def variant_class() -> pl.Expr:
    snv = (pl.col(EA).str.len_bytes() == 1) & (pl.col(NEA).str.len_bytes() == 1)
    palindromic = snv & (pl.col(EA).replace_strict(COMPLEMENT, default=None) == pl.col(NEA))
    return (
        pl.when(palindromic)
        .then(pl.lit("palindromic_snv"))
        .when(snv)
        .then(pl.lit("other_snv"))
        .otherwise(pl.lit("indel"))
        .alias("class")
    )


def read(path) -> pl.DataFrame:
    return pl.read_parquet(path, columns=COLUMNS).with_columns(
        pl.col(EA).cast(pl.String), pl.col(NEA).cast(pl.String)
    )


def report(record: Record) -> None:
    (name,) = record.file_md5
    old, new = read(record.old(name)), read(record.new(name))
    print(f"\n== {record.asset_id}: old={old.height:,} new={new.height:,}")
    for label, frame in (
        ("only in old (dropped)", old.join(new, on=GWASLAB_SNPID_COL, how="anti")),
        ("only in new (added)", new.join(old, on=GWASLAB_SNPID_COL, how="anti")),
    ):
        print(label)
        print(frame.group_by(variant_class()).len().sort("class"))
    shared = old.join(new, on=GWASLAB_SNPID_COL, suffix="_new")
    sign = shared.filter(
        pl.col(GWASLAB_BETA_COL).sign() != pl.col(GWASLAB_BETA_COL + "_new").sign()
    ).with_columns((pl.col(EA) != pl.col(EA + "_new")).alias("ea_changed"))
    print("shared rows whose BETA sign changed")
    print(sign.group_by(variant_class(), "ea_changed").len().sort("class", "ea_changed"))
    print("of those, a sample:")
    print(sign.head(5))


def main() -> None:
    with pl.Config(tbl_rows=20, tbl_cols=12, tbl_width_chars=200):
        for record in load_records():
            if record.task_type == "GenomeReferenceHarmonizationTask":
                report(record)


if __name__ == "__main__":
    main()
