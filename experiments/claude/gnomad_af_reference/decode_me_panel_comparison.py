"""
Validation step 3: harmonize DecodeME against four panel/ancestry choices and compare.

Choices: 1000 Genomes "eur", gnomAD v2.1.1 "nfe" and "nfe_nwe", Pan-UKBB "ukb_eur". For each,
the production trust decision and per-chromosome resolution run on DecodeME's
pre-harmonization table. Reported: the trust decision, drop-reason counts (kept, NOT_IN_PANEL,
AF_MISMATCH, ...), per-variant agreement of the outcome (kept orientation or drop reason)
between choices, and palindromes kept or flipped on a chosen-ancestry AF of exactly 0.

Usage (all four choices; needs the gnomAD v2.1.1 panel, a 14 h build):
    pixi r python -m experiments.claude.gnomad_af_reference.decode_me_panel_comparison \
        2>&1 | tee experiments/claude/gnomad_af_reference/decode_me_panel_comparison.log
Pass choice labels to compare a subset, building only the panels they need:
    pixi r python -m experiments.claude.gnomad_af_reference.decode_me_panel_comparison \
        1000g_eur pan_ukbb_eur \
        2>&1 | tee experiments/claude/gnomad_af_reference/decode_me_panel_comparison_1000g_pan_ukbb.log
"""

import itertools
import sys
from collections.abc import Sequence

import polars as pl
from attrs import frozen

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
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_allele_frequencies import (
    PAN_UKBB_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    PanelTable,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    resolve_chromosome_rows,
    resolve_panel_af_col,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    decide_trust,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.allele_frequency_panel_constants import PanelAncestry
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
OUTCOME_COL = "outcome"
KEPT = "kept"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
# The DecodeME rsID-assignment chain runs with palindromes dropped when unresolved.
OPTIONS = GenomeReferenceHarmonizationOptions()


@frozen(slots=True)
class PanelChoice:
    label: str
    panel_task: Task
    ancestry: PanelAncestry


CHOICES = [
    PanelChoice("1000g_eur", THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES, "eur"),
    PanelChoice("gnomad_nfe", GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES, "nfe"),
    PanelChoice(
        "gnomad_nfe_nwe", GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES, "nfe_nwe"
    ),
    PanelChoice("pan_ukbb_eur", PAN_UKBB_HG19_ALLELE_FREQUENCIES, "ukb_eur"),
]


def _file(asset: Asset) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def resolve_all(
    choice: PanelChoice, sumstats: pl.LazyFrame, fasta: IndexedFasta, panel: PanelTable
) -> pl.DataFrame:
    """Every input row, keyed by SNPID, with its outcome: oriented EA/NEA or a drop reason."""
    rules = resolve_column_rules(columns=sumstats.collect_schema().names(), extra=())
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
    evidence = count_trust_evidence_genome_wide(
        sumstats, chromosomes, fasta, panel, OPTIONS
    )
    trusted = decide_trust(evidence, OPTIONS)
    print(f"[{choice.label}] trusted={trusted} suspicious={evidence.suspicious}")
    parts = []
    for chrom in chromosomes:
        context = ChromosomeContext(
            chrom=chrom, fasta=fasta, trusted=trusted, rules=rules, options=OPTIONS
        )
        resolved = resolve_chromosome_rows(sumstats, context, panel)
        parts.append(
            resolved.select(
                GWASLAB_SNPID_COL,
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                EA,
                NEA,
                pl.coalesce(
                    pl.col(DROP_REASON_COL),
                    pl.lit(KEPT + ":") + pl.col(EA) + pl.lit("/") + pl.col(NEA),
                ).alias(OUTCOME_COL),
            )
        )
    return pl.concat(parts)


def report_drop_reasons(label: str, outcomes: pl.DataFrame) -> None:
    reason = (
        pl.when(pl.col(OUTCOME_COL).str.starts_with(KEPT))
        .then(pl.lit(KEPT))
        .otherwise(pl.col(OUTCOME_COL))
    )
    print(f"\n== [{label}] outcome counts")
    with pl.Config(tbl_rows=30):
        print(
            outcomes.group_by(reason.alias("reason")).len().sort("len", descending=True)
        )


def report_agreement(outcomes: dict[str, pl.DataFrame]) -> None:
    print("\n== pairwise outcome disagreement (rows whose outcome differs)")
    for left, right in itertools.combinations(outcomes, 2):
        joined = (
            outcomes[left]
            .select(GWASLAB_SNPID_COL, OUTCOME_COL)
            .join(
                outcomes[right].select(GWASLAB_SNPID_COL, OUTCOME_COL),
                on=GWASLAB_SNPID_COL,
                suffix="_right",
            )
        )
        differing = joined.filter(pl.col(OUTCOME_COL) != pl.col(OUTCOME_COL + "_right"))
        print(f"{left} vs {right}: {differing.height} of {joined.height}")
        both_kept = differing.filter(
            pl.col(OUTCOME_COL).str.starts_with(KEPT)
            & pl.col(OUTCOME_COL + "_right").str.starts_with(KEPT)
        )
        print(
            f"{left} vs {right}: kept by both in opposite orientations: {both_kept.height}"
        )
        with pl.Config(tbl_rows=10, fmt_str_lengths=40):
            print(
                differing.group_by(OUTCOME_COL, OUTCOME_COL + "_right")
                .len()
                .sort("len", descending=True)
                .head(10)
            )


def report_zero_af_palindromes(
    label: str, outcomes: pl.DataFrame, panel: PanelTable
) -> None:
    palindromic = (pl.col(EA).str.len_bytes() == 1) & (
        pl.col(EA).replace_strict(COMPLEMENT, default=None) == pl.col(NEA)
    )
    resolved = outcomes.filter(palindromic & pl.col(OUTCOME_COL).str.starts_with(KEPT))
    panel_rows = pl.scan_parquet(panel.path).select(
        GWASLAB_CHROM_COL,
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
        PANEL_REF_COL,
        PANEL_ALT_COL,
        pl.col(panel.af_col).alias("panel_af"),
    )
    matched = (
        resolved.lazy()
        .with_columns(
            pl.col(GWASLAB_CHROM_COL).cast(pl.Int32),
            pl.col(GWASLAB_POS_COL).cast(pl.Int64),
        )
        .join(panel_rows, on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL])
        .filter(
            (
                (pl.col(PANEL_REF_COL) == pl.col(NEA))
                & (pl.col(PANEL_ALT_COL) == pl.col(EA))
            )
            | (
                (pl.col(PANEL_REF_COL) == pl.col(EA))
                & (pl.col(PANEL_ALT_COL) == pl.col(NEA))
            )
        )
        .collect(engine="streaming")
    )
    print(
        f"\n== [{label}] resolved palindromes: {resolved.height}; "
        f"decided by a panel AF of exactly 0: {matched.filter(pl.col('panel_af') == 0).height}"
    )


def select_choices(labels: Sequence[str]) -> list[PanelChoice]:
    """The named choices, in CHOICES order; all of them when no label is given."""
    known = [choice.label for choice in CHOICES]
    unknown = sorted(set(labels) - set(known))
    assert not unknown, f"unknown choices {unknown}; known: {known}"
    return [choice for choice in CHOICES if not labels or choice.label in labels]


def main(labels: Sequence[str]) -> None:
    choices = select_choices(labels)
    pre_task = (
        DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.pre_harmonization_table_task
    )
    panel_tasks = list(
        {choice.panel_task.asset_id: choice.panel_task for choice in choices}.values()
    )
    assets = DEFAULT_RUNNER.run([pre_task, UCSC_HG19_INDEXED_FASTA, *panel_tasks])
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(
        assets[pre_task.asset_id], pre_task.meta, IdentityPipe()
    )
    assert GWASLAB_SNPID_COL in sumstats.collect_schema().names(), (
        f"the agreement join needs {GWASLAB_SNPID_COL} in the pre-harmonization table"
    )
    outcomes: dict[str, pl.DataFrame] = {}
    for choice in choices:
        panel = PanelTable(
            path=_file(assets[choice.panel_task.asset_id]).path,
            af_col=resolve_panel_af_col(choice.panel_task, choice.ancestry),
        )
        outcomes[choice.label] = resolve_all(choice, sumstats, fasta, panel)
        report_drop_reasons(choice.label, outcomes[choice.label])
        report_zero_af_palindromes(choice.label, outcomes[choice.label], panel)
    report_agreement(outcomes)


if __name__ == "__main__":
    main(sys.argv[1:])
