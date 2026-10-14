"""
How closely do DecodeME's resolved palindromic SNVs match their panel frequency?

Palindromes are resolved by gwaslab's side-of-0.5 rule: summary-statistic MAF and panel MAF
both at most 0.4, then kept if EAF and panel AF lie on the same side of 0.5 and strand-flipped
otherwise. Unlike the ambiguous-indel rule, there is no distance check, so a panel AF of 0
(or 0.01) settles a variant whose EAF is 0.38. This script measures, for each of the four
panel choices of decode_me_panel_comparison:

- the residual |EAF - panel AF| after harmonization (EAF in the output orientation, panel AF
  for REF=NEA, ALT=EA), binned;
- for palindromes decided by a panel AF of exactly 0 or 1, the distribution of their MAF;
- per residual bin, how often another choice that also resolved the variant chose the
  opposite strand. If a large residual flags unreliable decisions, disagreement rises with it.

Usage (needs all four panels):
    pixi r python -m experiments.claude.gnomad_af_reference.palindrome_eaf_distribution \
        2>&1 | tee experiments/claude/gnomad_af_reference/palindrome_eaf_distribution.log
"""

import itertools

import polars as pl

from experiments.claude.gnomad_af_reference.decode_me_panel_comparison import (
    COMPLEMENT,
    OPTIONS,
    PanelChoice,
    select_choices,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    PanelTable,
    ParquetPanelLoader,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    resolve_chromosome_rows,
    resolve_panel_af_col,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
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
EAF = GWASLAB_EFFECT_ALLELE_FREQ_COL
RESIDUAL_BREAKS = [0.02, 0.05, 0.1, 0.2, 0.3]
MAF_BREAKS = [0.001, 0.01, 0.05, 0.1, 0.2, 0.3]
FREQUENCY_TOLERANCE = 1e-9


def _file(asset: Asset) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def is_palindromic_snv() -> pl.Expr:
    return (pl.col(EA).str.len_bytes() == 1) & (
        pl.col(EA).replace_strict(COMPLEMENT, default=None) == pl.col(NEA)
    )


def resolved_palindromes(
    choice: PanelChoice, sumstats: pl.LazyFrame, fasta: IndexedFasta, panel: PanelTable
) -> pl.DataFrame:
    """Kept palindromic SNVs: input and output EAF, the strand decision and the panel AF."""
    rules = resolve_column_rules(columns=sumstats.collect_schema().names(), extra=())
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
    evidence = count_trust_evidence_genome_wide(
        sumstats, chromosomes, fasta, panel, OPTIONS
    )
    assert not decide_trust(evidence, OPTIONS), (
        f"{choice.label}: a trusted table resolves no palindromes by frequency"
    )
    parts = []
    for chrom in chromosomes:
        context = ChromosomeContext(
            chrom=chrom, fasta=fasta, trusted=False, rules=rules, options=OPTIONS
        )
        inputs = (
            sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
            .select(
                GWASLAB_SNPID_COL,
                pl.col(EA).cast(pl.String).str.to_uppercase().alias("ea_in"),
                pl.col(EAF).cast(pl.Float64).alias("eaf_in"),
            )
            .collect()
        )
        kept = (
            resolve_chromosome_rows(sumstats, context, panel)
            .with_columns(pl.col(EA).cast(pl.String), pl.col(NEA).cast(pl.String))
            .filter(pl.col(DROP_REASON_COL).is_null() & is_palindromic_snv())
        )
        records = ParquetPanelLoader(
            panel_path=panel.path, af_col=panel.af_col, chrom=chrom
        )(kept[GWASLAB_POS_COL])
        parts.append(
            kept.select(
                GWASLAB_SNPID_COL,
                GWASLAB_POS_COL,
                EA,
                NEA,
                pl.col(EAF).cast(pl.Float64).alias("eaf_out"),
            )
            .with_columns(pl.col(GWASLAB_POS_COL).cast(pl.Int64))
            .join(inputs, on=GWASLAB_SNPID_COL)
            .join(
                records.select(
                    GWASLAB_POS_COL,
                    pl.col(PANEL_REF_COL).alias(NEA),
                    pl.col(PANEL_ALT_COL).alias(EA),
                    pl.col(PANEL_AF_COL).cast(pl.Float64).alias("panel_af"),
                ),
                on=[GWASLAB_POS_COL, EA, NEA],
                how="left",
            )
        )
    frame = pl.concat(parts)
    assert frame["panel_af"].null_count() == 0, (
        f"{choice.label}: a resolved palindrome has no panel record"
    )
    # A swap moves EA and replaces EAF by 1 - EAF; a strand flip replaces EAF alone. So the
    # strand was flipped exactly when the EAF changed without the alleles swapping, or vice versa.
    eaf_changed = (pl.col("eaf_out") - pl.col("eaf_in")).abs() > FREQUENCY_TOLERANCE
    alleles_swapped = pl.col(EA) != pl.col("ea_in")
    return frame.select(
        GWASLAB_SNPID_COL,
        (eaf_changed ^ alleles_swapped).alias("strand_flip"),
        "eaf_out",
        "panel_af",
        (pl.col("eaf_out") - pl.col("panel_af")).abs().alias("residual"),
        pl.min_horizontal(pl.col("eaf_out"), 1 - pl.col("eaf_out")).alias("maf"),
    )


def binned(column: str, breaks: list[float]) -> pl.Expr:
    return pl.col(column).cut(breaks, left_closed=False).alias(column + "_bin")


def report_choice(label: str, frame: pl.DataFrame) -> None:
    print(f"\n== [{label}] resolved palindromes: {frame.height}")
    print(
        frame.select(
            pl.col("strand_flip").sum().alias("strand_flipped"),
            *[
                pl.col("residual").quantile(q).round(4).alias(f"residual_q{q}")
                for q in [0.5, 0.9, 0.99, 0.999]
            ],
            pl.col("residual").max().round(4).alias("residual_max"),
        )
    )
    print(f"[{label}] residual |EAF - panel AF| by bin")
    print(
        frame.group_by(binned("residual", RESIDUAL_BREAKS))
        .agg(pl.len(), pl.col("strand_flip").sum().alias("strand_flipped"))
        .sort("residual_bin")
    )
    monomorphic = frame.filter(
        (pl.col("panel_af") == 0) | (pl.col("panel_af") == 1)
    )
    print(
        f"[{label}] decided by a panel AF of exactly 0 or 1: {monomorphic.height}; "
        "their DecodeME MAF by bin"
    )
    print(
        monomorphic.group_by(binned("maf", MAF_BREAKS))
        .agg(pl.len(), pl.col("strand_flip").sum().alias("strand_flipped"))
        .sort("maf_bin")
    )


def report_cross_choice_disagreement(frames: dict[str, pl.DataFrame]) -> None:
    """Per residual bin of the first choice: how often a second choice flips the other way."""
    print("\n== strand disagreement between choices that both resolved a palindrome")
    for left, right in itertools.permutations(frames, 2):
        joined = frames[left].join(
            frames[right].select(GWASLAB_SNPID_COL, "strand_flip"),
            on=GWASLAB_SNPID_COL,
            suffix="_other",
        )
        print(f"\n[{left}] residual bin vs disagreement with [{right}]")
        print(
            joined.group_by(binned("residual", RESIDUAL_BREAKS))
            .agg(
                pl.len().alias("both_resolved"),
                (pl.col("strand_flip") != pl.col("strand_flip_other"))
                .sum()
                .alias("opposite_strand"),
            )
            .with_columns(
                (pl.col("opposite_strand") / pl.col("both_resolved"))
                .round(5)
                .alias("fraction")
            )
            .sort("residual_bin")
        )


def main() -> None:
    choices = select_choices([])
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
    frames: dict[str, pl.DataFrame] = {}
    with pl.Config(tbl_rows=20, tbl_cols=12, tbl_width_chars=200):
        for choice in choices:
            panel = PanelTable(
                path=_file(assets[choice.panel_task.asset_id]).path,
                af_col=resolve_panel_af_col(choice.panel_task, choice.ancestry),
            )
            frames[choice.label] = resolved_palindromes(choice, sumstats, fasta, panel)
            report_choice(choice.label, frames[choice.label])
        report_cross_choice_disagreement(frames)


if __name__ == "__main__":
    main()
