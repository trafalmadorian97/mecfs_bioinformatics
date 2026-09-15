"""
Write a markdown report of the 11 ambiguous indels that V3 resolves to "swap" on build-38
DecodeME at distance 0.02 and margin 0.3 (listed in inspect_v3_wrong_choices.log).

For each variant the report gives, as CHROM/POS/REF/ALT/ALT frequency:
- the raw regenie record (ALLELE0 as REF, ALLELE1 as ALT, A1FREQ);
- the build-38 gwaslab table (NEA as REF, EA as ALT, EAF);
- the 1000 Genomes EUR hg38 panel record for the source reading (REF = NEA, ALT = EA);
- the panel record for the opposite reading (REF = EA, ALT = NEA);
- every other panel record at that position, and the hg38 sequence around it.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.v3_wrong_choices_report \
    2>&1 | tee experiments/claude/genome_reference_harmonization/v3_wrong_choices_report.log
"""

from pathlib import Path

import polars as pl
import pysam
from attrs import frozen

from experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules import (
    BUILD_38_TABLE,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.filtered_snps_gwas_1 import (
    DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    FASTA_FILENAME,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_COL,
    REGENIE_ALLELE0_COL,
    REGENIE_ALLELE1_COL,
    REGENIE_CHROM_COL,
    REGENIE_GENPOS_COL,
    REGENIE_ID_COL,
)

REPORT_PATH = Path(
    "experiments/claude/genome_reference_harmonization/v3_wrong_choices.md"
)
FLANK_BASES = 20
ABSENT = "absent"


@frozen
class Site:
    chrom: int
    pos: int


# From inspect_v3_wrong_choices.log.
WRONG_SITES = [
    Site(2, 126390673),
    Site(2, 224984840),
    Site(2, 228597662),
    Site(11, 947010),
    Site(11, 19321212),
    Site(11, 47561137),
    Site(12, 57595702),
    Site(12, 64250003),
    Site(13, 32388121),
    Site(15, 69557108),
    Site(19, 43995406),
]


def _file(asset: Asset) -> Path:
    assert isinstance(asset, FileAsset)
    return asset.path


def _at(frame: pl.LazyFrame, chrom_col: str, pos_col: str, site: Site) -> pl.DataFrame:
    return frame.filter(
        (pl.col(chrom_col) == site.chrom) & (pl.col(pos_col) == site.pos)
    ).collect()


def _row(label: str, chrom: int, pos: int, ref: str, alt: str, freq: str) -> str:
    return f"| {label} | {chrom} | {pos} | {ref} | {alt} | {freq} |"


def _panel_af(panel: pl.DataFrame, ref: str, alt: str) -> str:
    match = panel.filter(
        (pl.col(PANEL_REF_COL) == ref) & (pl.col(PANEL_ALT_COL) == alt)
    )
    assert match.height <= 1
    return ABSENT if match.height == 0 else f"{match[PANEL_AF_COL][0]:.4f}"


def _section(
    site: Site,
    raw: pl.DataFrame,
    table: pl.DataFrame,
    panel: pl.DataFrame,
    genome: pysam.FastaFile,
) -> list[str]:
    assert table.height == 1, f"expected one build-38 row at {site}"
    ea = str(table[GWASLAB_EFFECT_ALLELE_COL][0])
    nea = str(table[GWASLAB_NON_EFFECT_ALLELE_COL][0])
    start = site.pos - 1
    before = genome.fetch(f"chr{site.chrom}", start - FLANK_BASES, start).upper()
    after = genome.fetch(f"chr{site.chrom}", start, start + FLANK_BASES).upper()
    lines = [
        f"## {site.chrom}:{site.pos}",
        "",
        f"hg38 sequence (position {site.pos} is the first base after the bar): "
        f"{before}|{after}",
        "",
        "| Source | CHROM | POS | REF | ALT | ALT frequency |",
        "|---|---|---|---|---|---|",
    ]
    for record in raw.iter_rows(named=True):
        lines.append(
            _row(
                f"DecodeME raw regenie (ID {record[REGENIE_ID_COL]}; ALLELE0, ALLELE1, A1FREQ)",
                site.chrom,
                site.pos,
                record[REGENIE_ALLELE0_COL],
                record[REGENIE_ALLELE1_COL],
                f"{record[REGENIE_A1FREQ_COL]:.4f}",
            )
        )
    lines.append(
        _row(
            f"DecodeME build-38 gwaslab table (SNPID {table[GWASLAB_SNPID_COL][0]}; NEA, EA, EAF)",
            site.chrom,
            site.pos,
            nea,
            ea,
            f"{table[GWASLAB_EFFECT_ALLELE_FREQ_COL][0]:.4f}",
        )
    )
    lines.append(
        _row(
            "1KG EUR panel, source reading (REF = NEA, ALT = EA)",
            site.chrom,
            site.pos,
            nea,
            ea,
            _panel_af(panel, ref=nea, alt=ea),
        )
    )
    lines.append(
        _row(
            "1KG EUR panel, opposite reading (REF = EA, ALT = NEA)",
            site.chrom,
            site.pos,
            ea,
            nea,
            _panel_af(panel, ref=ea, alt=nea),
        )
    )
    others = panel.filter(
        ~(
            (pl.col(PANEL_REF_COL) == nea) & (pl.col(PANEL_ALT_COL) == ea)
            | (pl.col(PANEL_REF_COL) == ea) & (pl.col(PANEL_ALT_COL) == nea)
        )
    )
    for record in others.iter_rows(named=True):
        lines.append(
            _row(
                "1KG EUR panel, other record at this position",
                site.chrom,
                site.pos,
                record[PANEL_REF_COL],
                record[PANEL_ALT_COL],
                f"{record[PANEL_AF_COL]:.4f}",
            )
        )
    lines.append("")
    return lines


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [
            DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
            BUILD_38_TABLE,
            UCSC_HG38_INDEXED_FASTA,
            THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
        ]
    )
    raw = pl.scan_parquet(_file(assets[DECODE_ME_FILTER_SNPS_GWAS_1_TASK.asset_id]))
    table = pl.scan_parquet(_file(assets[BUILD_38_TABLE.asset_id]))
    panel = pl.scan_parquet(
        _file(assets[THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES.asset_id])
    )
    fasta_asset = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    genome = pysam.FastaFile(str(fasta_asset.path / FASTA_FILENAME))
    lines = [
        "# V3 wrong choices: 11 ambiguous indels on build-38 DecodeME",
        "",
        "Ambiguous indels whose alleles both match hg38, resolved to swap by the stringent",
        "rules at distance 0.02 and margin 0.3, although build 38 is trusted (its source",
        "orientation is taken as truth). Frequencies are of the ALT allele. The panel is",
        "1000 Genomes 30x EUR on hg38, as built by ReferencePanelAlleleFrequencyTask.",
        "",
        "Generated by v3_wrong_choices_report.py.",
        "",
    ]
    for site in WRONG_SITES:
        lines.extend(
            _section(
                site,
                raw=_at(raw, REGENIE_CHROM_COL, REGENIE_GENPOS_COL, site),
                table=_at(table, GWASLAB_CHROM_COL, GWASLAB_POS_COL, site).with_columns(
                    pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String),
                    pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String),
                ),
                panel=_at(panel, GWASLAB_CHROM_COL, GWASLAB_POS_COL, site),
                genome=genome,
            )
        )
    REPORT_PATH.write_text("\n".join(lines))
    print(f"wrote {REPORT_PATH}")


if __name__ == "__main__":
    main()
