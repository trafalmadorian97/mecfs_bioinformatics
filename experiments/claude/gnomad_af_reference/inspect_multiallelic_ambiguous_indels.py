"""
DecodeME ambiguous indels at positions where gnomAD has several records.

An ambiguous indel is an input row whose EA and NEA both match the hg19 reference at POS
(for example NEA=A, EA=AT where the genome reads AT...). It has two readings: keep (REF=NEA,
ALT=EA: an insertion) and flip (REF=EA, ALT=NEA: a deletion). gnomAD splits multiallelic
sites into biallelic records, so both readings can exist as records for two different
variants at one position. The harmonizer then picks the reading whose AF is closer to EAF.

This script prints example positions with every gnomAD v2.1.1 record there (all groups'
AF, unfiltered), and counts DecodeME ambiguous indels by how many gnomAD records sit at
their position, which readings have a usable record (0 < AF_nfe < 1), and what the
ambiguous-indel rule decides with gnomAD nfe and with 1000 Genomes EUR.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.inspect_multiallelic_ambiguous_indels \
        2>&1 | tee experiments/claude/gnomad_af_reference/inspect_multiallelic_ambiguous_indels.log
"""

from pathlib import Path

import polars as pl
import pysam

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_BOTH,
    classify_alleles,
    prepare_alleles,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    INDEL_DROP_REASON_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    gwaslab_code_to_contig_name,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
POS = GWASLAB_POS_COL
OPTIONS = GenomeReferenceHarmonizationOptions()
GNOMAD_SHOWN_COLUMNS = ["AF_nfe", "AN_nfe", "AF_nfe_nwe", "AF_afr", "AF_eas", "AF_amr"]
EXAMPLES_PER_KIND = 4
_FLANK = 10


def _file(asset: Asset) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def ambiguous_indels(
    sumstats: pl.LazyFrame, fasta: IndexedFasta, chrom: int
) -> pl.DataFrame:
    """DecodeME rows of one chromosome whose EA and NEA both match the genome."""
    rows = prepare_alleles(
        sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
        .select(GWASLAB_SNPID_COL, POS, EA, NEA, GWASLAB_EFFECT_ALLELE_FREQ_COL)
        .collect(engine="streaming")
    ).filter(valid_alleles_expr())
    classified = classify_alleles(
        rows, fasta=fasta, chrom=chrom, max_gather_bytes=OPTIONS.max_gather_bytes
    )
    return classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH).select(
        GWASLAB_SNPID_COL,
        POS,
        EA,
        NEA,
        pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64),
    )


def panel_records(
    path: Path, chrom: int, positions: pl.Series, columns: list[str]
) -> pl.DataFrame:
    """Every panel record at the positions, unfiltered."""
    return (
        pl.scan_parquet(path)
        .filter(pl.col(GWASLAB_CHROM_COL) == chrom)
        .filter(pl.col(POS).is_in(positions.unique().implode()))
        .select(pl.col(POS).cast(pl.Int64), PANEL_REF_COL, PANEL_ALT_COL, *columns)
        .collect()
    )


def rule_decision(rows: pl.DataFrame, records: pl.DataFrame, af_col: str) -> pl.Series:
    """What the ambiguous-indel rule decides for each row against these records."""
    panel = records.select(
        POS, PANEL_REF_COL, PANEL_ALT_COL, pl.col(af_col).alias(PANEL_AF_COL)
    ).filter(pl.col(PANEL_AF_COL).is_not_null())
    decided = decide_ambiguous_indels(
        rows.select(POS, EA, NEA, GWASLAB_EFFECT_ALLELE_FREQ_COL), panel, OPTIONS
    )
    return decided.select(
        pl.coalesce(INDEL_ACTION_COL, INDEL_DROP_REASON_COL)
    ).to_series()


def reading_flags(rows: pl.DataFrame, records: pl.DataFrame) -> pl.DataFrame:
    """Records at the position, and whether each reading has a usable gnomAD nfe record."""
    usable = records.filter((pl.col("AF_nfe") > 0) & (pl.col("AF_nfe") < 1))
    keep = usable.select(
        POS,
        pl.col(PANEL_REF_COL).alias(NEA),
        pl.col(PANEL_ALT_COL).alias(EA),
        pl.lit(True).alias("has_keep"),
    )
    flip = usable.select(
        POS,
        pl.col(PANEL_REF_COL).alias(EA),
        pl.col(PANEL_ALT_COL).alias(NEA),
        pl.lit(True).alias("has_flip"),
    )
    counts = records.group_by(POS).len("gnomad_records_at_pos")
    mirrored = rows.select(
        POS, pl.col(NEA).alias(EA), pl.col(EA).alias(NEA), pl.lit(True).alias("mirror")
    )
    return (
        rows.join(counts, on=POS, how="left")
        .join(keep, on=[POS, EA, NEA], how="left")
        .join(flip, on=[POS, EA, NEA], how="left")
        .join(mirrored, on=[POS, EA, NEA], how="left")
        .with_columns(
            pl.col("gnomad_records_at_pos").fill_null(0),
            pl.col("has_keep").fill_null(False),
            pl.col("has_flip").fill_null(False),
            pl.col("mirror").fill_null(False).alias("decodeme_has_mirrored_row"),
        )
        .drop("mirror")
    )


def print_example(
    row: dict, chrom: int, gnomad: pl.DataFrame, kg: pl.DataFrame, genome: pysam.FastaFile
) -> None:
    pos = row[POS]
    contig = "chr" + gwaslab_code_to_contig_name(chrom)
    context = genome.fetch(contig, pos - 1 - _FLANK, pos - 1 + _FLANK + 1)
    print(
        f"\n### {row[GWASLAB_SNPID_COL]} {contig}:{pos}  DecodeME EA={row[EA]} "
        f"NEA={row[NEA]} EAF={row[GWASLAB_EFFECT_ALLELE_FREQ_COL]:.4f}  "
        f"rule: gnomad_nfe={row['gnomad_nfe_decision']} 1000g={row['1000g_decision']}"
    )
    print(
        f"hg19 {pos - _FLANK}-{pos + _FLANK}: "
        f"{context[:_FLANK]}[{context[_FLANK]}]{context[_FLANK + 1 :]}"
    )
    print("gnomAD v2.1.1 records at the position:")
    print(gnomad.filter(pl.col(POS) == pos))
    print("1000 Genomes EUR records at the position:")
    print(kg.filter(pl.col(POS) == pos))


def main() -> None:
    pre_task = (
        DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.pre_harmonization_table_task
    )
    assets = DEFAULT_RUNNER.run(
        [
            pre_task,
            UCSC_HG19_INDEXED_FASTA,
            GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
            THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
        ]
    )
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    genome = pysam.FastaFile(str(fasta.fasta_path))
    gnomad_path = _file(
        assets[GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES.asset_id]
    ).path
    kg_path = _file(
        assets[THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES.asset_id]
    ).path
    sumstats = scan_sumstats_as_polars(
        assets[pre_task.asset_id], pre_task.meta, IdentityPipe()
    )
    chromosomes = (
        sumstats.select(pl.col(GWASLAB_CHROM_COL).unique().sort())
        .collect()
        .to_series()
        .to_list()
    )
    parts = []
    examples: dict[str, int] = {}
    with pl.Config(tbl_rows=20, tbl_cols=12, tbl_width_chars=200, fmt_str_lengths=40):
        for chrom in chromosomes:
            rows = ambiguous_indels(sumstats, fasta, chrom)
            gnomad = panel_records(gnomad_path, chrom, rows[POS], GNOMAD_SHOWN_COLUMNS)
            kg = panel_records(kg_path, chrom, rows[POS], ["AF"])
            flagged = reading_flags(rows, gnomad).with_columns(
                rule_decision(rows, gnomad, "AF_nfe").alias("gnomad_nfe_decision"),
                rule_decision(rows, kg, "AF").alias("1000g_decision"),
            )
            parts.append(flagged)
            # Examples from the first chromosomes that have them, a few of each kind.
            for kind, condition in [
                ("both readings in gnomAD", pl.col("has_keep") & pl.col("has_flip")),
                (
                    "one reading, other records at the position",
                    (pl.col("has_keep") ^ pl.col("has_flip"))
                    & (pl.col("gnomad_records_at_pos") > 1),
                ),
            ]:
                for row in flagged.filter(condition).head(2).iter_rows(named=True):
                    if examples.get(kind, 0) >= EXAMPLES_PER_KIND:
                        break
                    examples[kind] = examples.get(kind, 0) + 1
                    print(f"\n======== example: {kind}")
                    print_example(row, chrom, gnomad, kg, genome)
    flagged = pl.concat(parts)
    multiallelic = pl.col("gnomad_records_at_pos") > 1
    with pl.Config(tbl_rows=60, tbl_width_chars=200, fmt_str_lengths=40):
        print(f"\n== DecodeME ambiguous indels: {flagged.height}")
        print(
            "\n== by gnomAD records at the position (any AF, any filter that reached the panel)"
        )
        print(
            flagged.group_by(
                pl.col("gnomad_records_at_pos").clip(upper_bound=4).alias("records (4=4+)")
            )
            .len()
            .sort("records (4=4+)")
        )
        print("\n== multiallelic (>1 gnomAD record) x usable readings x gnomAD nfe rule")
        print(
            flagged.group_by(
                multiallelic.alias("multiallelic"),
                "has_keep",
                "has_flip",
                "gnomad_nfe_decision",
            )
            .len()
            .sort("len", descending=True)
        )
        print("\n== multiallelic x 1000 Genomes EUR rule")
        print(
            flagged.group_by(multiallelic.alias("multiallelic"), "1000g_decision")
            .len()
            .sort("len", descending=True)
        )
        print("\n== gnomAD nfe rule x 1000 Genomes EUR rule, multiallelic sites only")
        print(
            flagged.filter(multiallelic)
            .group_by("gnomad_nfe_decision", "1000g_decision")
            .len()
            .sort("len", descending=True)
        )
        print(
            "\n== DecodeME also has the mirrored row (EA and NEA swapped) at the position"
        )
        print(
            flagged.group_by(
                multiallelic.alias("multiallelic"), "decodeme_has_mirrored_row"
            )
            .len()
            .sort("multiallelic", "decodeme_has_mirrored_row")
        )


if __name__ == "__main__":
    main()
