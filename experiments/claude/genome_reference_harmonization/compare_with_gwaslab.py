"""
V1/V2: compare genome-reference harmonization with gwaslab harmonization, row by row.

For DecodeME (liftover to 37) and Liu et al. 2023 IBD:
- build the pre-harmonization table and run GenomeReferenceHarmonizationTask on it;
- recompute per-chromosome drop reasons with the library;
- join the output with the cached gwaslab-harmonized table on SNPID;
- bucket every row and break the buckets down by variant type.

The cached gwaslab tables may predate the installed gwaslab version, because the build
cache does not key on code.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.compare_with_gwaslab \
    2>&1 | tee experiments/claude/genome_reference_harmonization/compare_with_gwaslab.log
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_dump_to_parquet import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats import (
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import (
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    reverse_complement_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    resolve_chromosome_rows,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    decide_trust,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
    GWASLAB_STATUS_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_CASES_COL,
    REGENIE_A1FREQ_CONTROLS_COL,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
BETA = GWASLAB_BETA_COL
SUFFIX = "_gwaslab"
STAT_TOLERANCE = 1e-6
BUCKET_COL = "bucket"
VARIANT_TYPE_COL = "variant_type"
INDEL_INFERENCE_FLIP_COL = "gwaslab_indel_inference_flip"
TYPE_INDEL = "indel"
TYPE_PALINDROMIC_SNV = "palindromic_snv"
TYPE_SNV = "snv"
TYPE_MNP = "mnp"
BUCKET_NEW_ONLY = "new_only"
BUCKET_GWASLAB_ONLY = "gwaslab_only"
BUCKET_IDENTICAL = "identical"
BUCKET_SAME_ORIENTATION_BETA_DIFFERS = "same_orientation_beta_differs"
BUCKET_OPPOSITE_ORIENTATION = "opposite_orientation"
BUCKET_OTHER = "other"
# gwaslab STATUS digit 7 value 4: stats flipped by indel inference (the bug signature)
GWASLAB_STATUS_DIGIT7_INDEL_FLIPPED = "4"
OPTIONS = GenomeReferenceHarmonizationOptions()


@frozen
class ComparisonCase:
    label: str
    pre_harmonization: Task
    gwaslab_harmonized: Task


@frozen
class DropReport:
    trusted: bool
    dropped: pl.DataFrame


CASES = [
    ComparisonCase(
        label="decode_me_gwas_1",
        pre_harmonization=GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
            asset_id="experiment_decode_me_gwas_1_pre_harmonization_table",
            sub_dir="processed",
        ),
        gwaslab_harmonized=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.dump_parquet_task,
    ),
    ComparisonCase(
        label="liu_et_al_2023_ibd",
        pre_harmonization=GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
            asset_id="experiment_liu_et_al_2023_ibd_pre_harmonization_table",
            sub_dir="processed",
        ),
        gwaslab_harmonized=LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
    ),
]


def _file(asset: Asset) -> Path:
    assert isinstance(asset, FileAsset)
    return asset.path


def _variant_type(ea: str, nea: str) -> pl.Expr:
    ea_length = pl.col(ea).str.len_bytes()
    nea_length = pl.col(nea).str.len_bytes()
    palindromic = (
        (ea_length == 1)
        & (nea_length == 1)
        & (reverse_complement_expr(nea) == pl.col(ea))
    )
    return (
        pl.when(ea_length != nea_length)
        .then(pl.lit(TYPE_INDEL))
        .when(palindromic)
        .then(pl.lit(TYPE_PALINDROMIC_SNV))
        .when(ea_length == 1)
        .then(pl.lit(TYPE_SNV))
        .otherwise(pl.lit(TYPE_MNP))
    )


def drop_report(case: ComparisonCase, assets: dict[str, Asset]) -> DropReport:
    pre_asset = assets[case.pre_harmonization.asset_id]
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    panel_path = _file(
        assets[THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES.asset_id]
    )
    sumstats = scan_sumstats_as_polars(
        pre_asset, case.pre_harmonization.meta, IdentityPipe()
    )
    fasta = IndexedFasta.open(fasta_asset.path)
    rules = resolve_column_rules(columns=sumstats.collect_schema().names(), extra=())
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
    counts = count_trust_evidence_genome_wide(sumstats, chromosomes, fasta, OPTIONS)
    trusted = decide_trust(counts, OPTIONS)
    print(f"[{case.label}] trust counts {counts} -> trusted={trusted}")
    parts = []
    for chrom in chromosomes:
        context = ChromosomeContext(
            chrom=chrom, fasta=fasta, trusted=trusted, rules=rules, options=OPTIONS
        )
        resolved = resolve_chromosome_rows(sumstats, context, panel_path)
        parts.append(
            resolved.filter(pl.col(DROP_REASON_COL).is_not_null()).select(
                GWASLAB_SNPID_COL,
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                EA,
                NEA,
                DROP_REASON_COL,
            )
        )
    return DropReport(trusted=trusted, dropped=pl.concat(parts))


def compare(
    case: ComparisonCase, new: pl.DataFrame, old: pl.DataFrame, drops: DropReport
) -> None:
    for name, frame in [("new", new), ("gwaslab", old)]:
        assert frame[GWASLAB_SNPID_COL].is_unique().all(), (
            f"{name} SNPIDs are not unique"
        )
    old = old.with_columns(pl.col(EA).cast(pl.String), pl.col(NEA).cast(pl.String))
    joined = new.join(
        old, on=GWASLAB_SNPID_COL, how="full", suffix=SUFFIX, coalesce=True
    )
    same_orientation = (pl.col(EA) == pl.col(EA + SUFFIX)) & (
        pl.col(NEA) == pl.col(NEA + SUFFIX)
    )
    opposite_orientation = (pl.col(EA) == pl.col(NEA + SUFFIX)) & (
        pl.col(NEA) == pl.col(EA + SUFFIX)
    )
    bucket = (
        pl.when(pl.col(EA + SUFFIX).is_null())
        .then(pl.lit(BUCKET_NEW_ONLY))
        .when(pl.col(EA).is_null())
        .then(pl.lit(BUCKET_GWASLAB_ONLY))
        .when(
            same_orientation
            & ((pl.col(BETA) - pl.col(BETA + SUFFIX)).abs() < STAT_TOLERANCE)
        )
        .then(pl.lit(BUCKET_IDENTICAL))
        .when(same_orientation)
        .then(pl.lit(BUCKET_SAME_ORIENTATION_BETA_DIFFERS))
        .when(
            opposite_orientation
            & ((pl.col(BETA) + pl.col(BETA + SUFFIX)).abs() < STAT_TOLERANCE)
        )
        .then(pl.lit(BUCKET_OPPOSITE_ORIENTATION))
        .otherwise(pl.lit(BUCKET_OTHER))
    )
    joined = joined.with_columns(
        bucket.alias(BUCKET_COL),
        pl.coalesce(
            _variant_type(EA, NEA), _variant_type(EA + SUFFIX, NEA + SUFFIX)
        ).alias(VARIANT_TYPE_COL),
        (
            pl.col(GWASLAB_STATUS_COL + SUFFIX).cast(pl.String).str.slice(6, 1)
            == GWASLAB_STATUS_DIGIT7_INDEL_FLIPPED
        ).alias(INDEL_INFERENCE_FLIP_COL),
    )
    with pl.Config(tbl_rows=60, tbl_cols=-1, tbl_width_chars=250):
        print(f"\n=== {case.label}: buckets by variant type ===")
        print(
            joined.group_by(BUCKET_COL, VARIANT_TYPE_COL)
            .len()
            .sort(BUCKET_COL, VARIANT_TYPE_COL)
        )
        print(
            f"\n=== {case.label}: opposite-orientation rows by gwaslab indel-inference flip ==="
        )
        print(
            joined.filter(pl.col(BUCKET_COL) == BUCKET_OPPOSITE_ORIENTATION)
            .group_by(VARIANT_TYPE_COL, INDEL_INFERENCE_FLIP_COL)
            .len()
        )
        gwaslab_only = joined.filter(pl.col(BUCKET_COL) == BUCKET_GWASLAB_ONLY).join(
            drops.dropped.select(GWASLAB_SNPID_COL, DROP_REASON_COL),
            on=GWASLAB_SNPID_COL,
            how="left",
        )
        print(f"\n=== {case.label}: rows only gwaslab kept, by our drop reason ===")
        print(
            gwaslab_only.group_by(DROP_REASON_COL, VARIANT_TYPE_COL)
            .len()
            .sort("len", descending=True)
        )
        print(f"\n=== {case.label}: rows only we kept ===")
        print(
            joined.filter(pl.col(BUCKET_COL) == BUCKET_NEW_ONLY)
            .group_by(VARIANT_TYPE_COL)
            .len()
        )
        for name in [
            BUCKET_SAME_ORIENTATION_BETA_DIFFERS,
            BUCKET_OTHER,
            BUCKET_NEW_ONLY,
            BUCKET_GWASLAB_ONLY,
        ]:
            print(f"\n--- {case.label}: examples of {name} ---")
            print(joined.filter(pl.col(BUCKET_COL) == name).head(8))
        for column in [REGENIE_A1FREQ_CASES_COL, REGENIE_A1FREQ_CONTROLS_COL]:
            if column in new.columns:
                differs = joined.filter(
                    (pl.col(BUCKET_COL) == BUCKET_IDENTICAL)
                    & (
                        (pl.col(column) - pl.col(column + SUFFIX)).abs()
                        > STAT_TOLERANCE
                    )
                )
                print(
                    f"\n{case.label}: identical rows whose {column} differs (gwaslab never flips it): {differs.height:,}"
                )


def main() -> None:
    for case in CASES:
        new_task = GenomeReferenceHarmonizationTask.create(
            asset_id=f"experiment_{case.label}_genome_reference_harmonized",
            sumstats_task=case.pre_harmonization,
            fasta_task=UCSC_HG19_INDEXED_FASTA,
            panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
            options=OPTIONS,
        )
        assets = dict(
            DEFAULT_RUNNER.run(
                [
                    new_task,
                    case.pre_harmonization,
                    case.gwaslab_harmonized,
                    UCSC_HG19_INDEXED_FASTA,
                    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
                ]
            )
        )
        new = pl.read_parquet(_file(assets[new_task.asset_id]))
        old = pl.read_parquet(_file(assets[case.gwaslab_harmonized.asset_id]))
        compare(case, new, old, drop_report(case, assets))


if __name__ == "__main__":
    main()
