"""
Context of the "*" and X:X rows that sit alone at their position in the deCODE reference file
(AP4M1): each with every row within 12 bp on either side and the GRCh38 reference bases.

  pixi r python experiments/claude/decode_proteomics_survey/explore_lone_shadow_rows.py
"""

import polars as pl
import pysam

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    FASTA_FILENAME,
)

FLANK = 12
N_EXAMPLES = 6


def read_rows(path) -> pl.DataFrame:  # noqa: ANN001
    ea, oa = pl.col("effectAllele"), pl.col("otherAllele")
    return pl.read_csv(
        path,
        separator="\t",
        columns=["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "Beta", "SE", "ImpMAF"],
        null_values=["."],
        quote_char=None,
        schema_overrides={"Pos": pl.Int64},
    ).with_columns(
        pl.when((ea == "*") | (oa == "*"))
        .then(pl.lit("star"))
        .when(ea == oa)
        .then(pl.lit("same_allele"))
        .when((ea.str.len_bytes() == 1) & (oa.str.len_bytes() == 1))
        .then(pl.lit("snv"))
        .otherwise(pl.lit("indel"))
        .alias("kind"),
        pl.len().over("Chrom", "Pos").alias("n_at_pos"),
    )


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, UCSC_HG38_INDEXED_FASTA])
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    fasta_dir = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(fasta_dir, DirectoryAsset)
    fasta = pysam.FastaFile(str(fasta_dir.path / FASTA_FILENAME))
    rows = read_rows(raw.path).filter(pl.col("Chrom") == "chr1")
    pl.Config.set_tbl_rows(40)
    pl.Config.set_tbl_width_chars(200)
    for kind in ("star", "same_allele"):
        lone = rows.filter((pl.col("kind") == kind) & (pl.col("n_at_pos") == 1))
        print(f"\n##### lone {kind} rows on chr1: {lone.height}")
        for pos in lone["Pos"].gather_every(max(lone.height // N_EXAMPLES, 1))[:N_EXAMPLES]:
            ref = fasta.fetch("chr1", pos - 1 - FLANK, pos + FLANK).upper()
            print(f"\n--- chr1:{pos}  reference {pos - FLANK}..{pos + FLANK}: "
                  f"{ref[:FLANK]} [{ref[FLANK]}] {ref[FLANK + 1:]}")
            print(
                rows.filter(pl.col("Pos").is_between(pos - FLANK, pos + FLANK)).select(
                    "Pos", "effectAllele", "otherAllele", "kind", "Beta", "SE", "ImpMAF"
                )
            )


if __name__ == "__main__":
    main()
