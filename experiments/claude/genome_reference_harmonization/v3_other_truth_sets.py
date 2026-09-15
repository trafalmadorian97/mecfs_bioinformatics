"""
Does the V3 ambiguous-indel failure seen on build-38 DecodeME recur in other fully trusted
summary statistics?

For one gwaslab Sumstats pickle in its source build (no liftover):
- pick the genome build (hg19 or hg38) whose FASTA matches more chromosome-1 SNVs;
- measure trust genome-wide with the default options;
- if trusted and EAF is present, resolve every ambiguous indel with the stringent rules
  over a few settings and rule variants, counting "swap" decisions as wrong;
- print the wrong choices at the strictest setting with all panel records at the position.

A calibration line compares EAF with panel AF for chromosome-1 non-palindromic SNVs
whose NEA is the reference base, to check that EAF is the effect-allele frequency in a
population like the panel's.

Run all candidates (one process each) with:
  for p in <pickles>; do
    pixi r python -m experiments.claude.genome_reference_harmonization.v3_other_truth_sets $p
  done 2>&1 | tee experiments/claude/genome_reference_harmonization/v3_other_truth_sets.log
"""

import sys
import tempfile
from pathlib import Path

import attrs
import gwaslab
import polars as pl
from attrs import frozen

from experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules import (
    BASE_OPTIONS,
    ChromosomeAmbiguousIndels,
    _ambiguous_indels,
)
from experiments.claude.genome_reference_harmonization.v3_rule_variants import (
    _both_records_mask,
    _drop_monomorphic,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_NEA_REF,
    IS_PALINDROMIC_SNV_COL,
    classify_alleles,
    prepare_alleles,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_KEEP,
    ACTION_SWAP,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    TrustCounts,
    decide_trust,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

KEY_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]
SETTINGS = ((0.02, 0.05), (0.02, 0.3), (0.1, 0.2))
STRICT_DISTANCE, STRICT_MARGIN = 0.02, 0.3
PANEL_RECORDS_COL = "panel_records_at_position"
MAX_EXAMPLES = 30
BUILD_PROBE_CHROMOSOME = 1


@frozen
class Build:
    name: str
    fasta_task: Task
    panel_task: Task


BUILDS = (
    Build(
        "hg19",
        UCSC_HG19_INDEXED_FASTA,
        THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
    ),
    Build(
        "hg38",
        UCSC_HG38_INDEXED_FASTA,
        THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
    ),
)


def _load_table(pickle_path: Path, parquet_path: Path) -> bool:
    data = gwaslab.load_pickle(str(pickle_path)).data
    columns = [
        c for c in [*KEY_COLUMNS, GWASLAB_EFFECT_ALLELE_FREQ_COL] if c in data.columns
    ]
    frame = pl.from_pandas(data[columns]).with_columns(
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String),
        pl.col(GWASLAB_CHROM_COL).cast(pl.Int64),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
    )
    del data
    n_rows = frame.height
    frame = frame.drop_nulls(KEY_COLUMNS)
    print(f"rows: {n_rows:,}; rows with null keys dropped: {n_rows - frame.height:,}")
    frame.write_parquet(parquet_path)
    return GWASLAB_EFFECT_ALLELE_FREQ_COL in columns


def _consistency(counts: TrustCounts) -> str:
    snvs = counts.consistent_snvs + counts.inconsistent_snvs
    indels = counts.consistent_indels + counts.inconsistent_indels
    return (
        f"SNV consistency {counts.consistent_snvs / max(snvs, 1):.5%} of {snvs:,}; "
        f"indel consistency {counts.consistent_indels / max(indels, 1):.5%} of {indels:,}"
    )


@frozen
class OpenedBuild:
    fasta: IndexedFasta
    panel_path: Path


def _open(build: Build) -> OpenedBuild:
    assets = DEFAULT_RUNNER.run([build.fasta_task, build.panel_task])
    fasta_asset = assets[build.fasta_task.asset_id]
    panel_asset = assets[build.panel_task.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset) and isinstance(
        panel_asset, FileAsset
    )
    return OpenedBuild(
        fasta=IndexedFasta.open(fasta_asset.path), panel_path=panel_asset.path
    )


def _print_eaf_calibration(
    sumstats: pl.LazyFrame, fasta: IndexedFasta, panel_path: Path
) -> None:
    chrom = BUILD_PROBE_CHROMOSOME
    rows = (
        sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
        .select(
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
            GWASLAB_EFFECT_ALLELE_FREQ_COL,
        )
        .collect()
    )
    valid = prepare_alleles(rows).filter(valid_alleles_expr())
    snvs = valid.filter(
        (pl.col(GWASLAB_EFFECT_ALLELE_COL).str.len_bytes() == 1)
        & (pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).str.len_bytes() == 1)
    )
    classified = classify_alleles(
        snvs, fasta=fasta, chrom=chrom, max_gather_bytes=BASE_OPTIONS.max_gather_bytes
    )
    checkable = classified.filter(
        (pl.col(ALLELE_CLASS_COL) == CLASS_NEA_REF) & ~pl.col(IS_PALINDROMIC_SNV_COL)
    )
    panel = (
        pl.scan_parquet(panel_path)
        .filter(pl.col(GWASLAB_CHROM_COL) == chrom)
        .select(
            pl.col(GWASLAB_POS_COL).cast(pl.Int64),
            pl.col(PANEL_REF_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
            pl.col(PANEL_ALT_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
            pl.col(PANEL_AF_COL).cast(pl.Float64),
        )
        .collect()
    )
    joined = checkable.join(
        panel,
        on=[GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL],
    )
    eaf = pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    af = pl.col(PANEL_AF_COL)
    summary = joined.select(
        pl.len().alias("n"),
        (eaf - af).abs().median().alias("median_abs_eaf_minus_af"),
        ((eaf - af).abs() <= 0.05).mean().alias("frac_within_0_05"),
        ((eaf - (1 - af)).abs() < (eaf - af).abs())
        .mean()
        .alias("frac_closer_to_1_minus_af"),
    ).row(0, named=True)
    print(f"EAF calibration (chr{chrom} non-palindromic SNVs in panel): {summary}")


def _wrong_examples(item: ChromosomeAmbiguousIndels, chrom: int) -> pl.DataFrame:
    options = attrs.evolve(
        BASE_OPTIONS,
        indel_max_af_distance=STRICT_DISTANCE,
        indel_min_af_margin=STRICT_MARGIN,
    )
    decisions = decide_ambiguous_indels(item.rows, item.panel, options)
    wrong = item.rows.filter(decisions[INDEL_ACTION_COL] == ACTION_SWAP)
    records = (
        item.panel.join(wrong.select(GWASLAB_POS_COL).unique(), on=GWASLAB_POS_COL)
        .group_by(GWASLAB_POS_COL)
        .agg(
            pl.format(
                "{}>{}:{}",
                pl.col(PANEL_REF_COL),
                pl.col(PANEL_ALT_COL),
                pl.col(PANEL_AF_COL).round(3),
            )
            .str.join("; ")
            .alias(PANEL_RECORDS_COL)
        )
    )
    return wrong.with_columns(pl.lit(chrom).alias(GWASLAB_CHROM_COL)).join(
        records, on=GWASLAB_POS_COL, how="left"
    )


def main() -> None:
    pickle_path = Path(sys.argv[1])
    print(f"\n######## {pickle_path.name}")
    with tempfile.TemporaryDirectory() as scratch:
        parquet_path = Path(scratch) / "table.parquet"
        has_eaf = _load_table(pickle_path, parquet_path)
        if not has_eaf:
            print("no EAF column: cannot test the stringent rules")
            return
        sumstats = pl.scan_parquet(parquet_path)
        consistent_by_build: dict[str, int] = {}
        for candidate in BUILDS:
            counts = count_trust_evidence_genome_wide(
                sumstats, [BUILD_PROBE_CHROMOSOME], _open(candidate).fasta, BASE_OPTIONS
            )
            consistent_by_build[candidate.name] = counts.consistent_snvs
            print(
                f"chr{BUILD_PROBE_CHROMOSOME} vs {candidate.name}: {_consistency(counts)}"
            )
        build = max(BUILDS, key=lambda candidate: consistent_by_build[candidate.name])
        opened = _open(build)
        fasta, panel_path = opened.fasta, opened.panel_path
        sumstats = sumstats.filter(pl.col(GWASLAB_CHROM_COL).is_in(list(fasta.entries)))
        chromosomes = chromosomes_to_harmonize(sumstats, fasta, BASE_OPTIONS)
        counts = count_trust_evidence_genome_wide(
            sumstats, chromosomes, fasta, BASE_OPTIONS
        )
        trusted = decide_trust(counts, BASE_OPTIONS)
        print(f"build {build.name}: {counts}")
        print(f"{_consistency(counts)} -> trusted={trusted}")
        if not trusted:
            return
        _print_eaf_calibration(sumstats, fasta, panel_path)
        per_chromosome = {
            chrom: _ambiguous_indels(sumstats, chrom, fasta, panel_path)
            for chrom in chromosomes
        }
        print(
            f"ambiguous indels: {sum(item.rows.height for item in per_chromosome.values()):,}"
        )
        variants = {
            "baseline": list(per_chromosome.values()),
            "drop_monomorphic": [
                attrs.evolve(item, panel=_drop_monomorphic(item.panel))
                for item in per_chromosome.values()
            ],
            "require_both_records": [
                ChromosomeAmbiguousIndels(
                    rows=item.rows.filter(_both_records_mask(item)), panel=item.panel
                )
                for item in per_chromosome.values()
            ],
        }
        for name, items in variants.items():
            for distance, margin in SETTINGS:
                options = attrs.evolve(
                    BASE_OPTIONS,
                    indel_max_af_distance=distance,
                    indel_min_af_margin=margin,
                )
                decisions = pl.concat(
                    [
                        decide_ambiguous_indels(item.rows, item.panel, options)
                        for item in items
                    ]
                )
                kept = int((decisions[INDEL_ACTION_COL] == ACTION_KEEP).sum())
                wrong = int((decisions[INDEL_ACTION_COL] == ACTION_SWAP).sum())
                print(
                    f"{name:22s} distance={distance} margin={margin}: kept={kept:,} wrong={wrong:,}"
                )
        examples = pl.concat(
            [_wrong_examples(item, chrom) for chrom, item in per_chromosome.items()]
        ).select(
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
            GWASLAB_EFFECT_ALLELE_FREQ_COL,
            PANEL_RECORDS_COL,
        )
        with pl.Config(
            tbl_rows=MAX_EXAMPLES, tbl_cols=-1, tbl_width_chars=250, fmt_str_lengths=120
        ):
            print(
                f"wrong choices (baseline, distance {STRICT_DISTANCE}, margin {STRICT_MARGIN}):"
            )
            print(examples.head(MAX_EXAMPLES))


if __name__ == "__main__":
    main()
