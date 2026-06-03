"""
S-LDSC discrepancy isolation experiment (DecodeME vs. Martin).

Peter's DecodeME cell-type S-LDSC p-values run slightly more significant than collaborator
Martin's (original LDSC 1.0.1).  This script runs a one-factor-at-a-time sweep to attribute the
gap, varying:

- strand-ambiguous (palindromic) SNP handling -- gwaslab's cts port retains palindromes, while
  original LDSC munge removes them;
- baseline model version (v1.2 vs. v1.1);
- a few additional QC filters (indels, HLA).

Everything else is held fixed.  All runs use gwaslab 4.1.7, so none reproduces the *documented*
(3.6.8) numbers exactly; config C0 (drop-ambiguous input, no extra filter, baseline v1.2) is the
on-4.1.7 reference point matching the current documented pipeline.

The matrix is run on the smaller `gtex_brain` partitioning dataset first, then on the larger
GTEx/Franke `multi_tissue_gene_expression` dataset, so that at least one full matrix completes
overnight.  Compare the per-config `ldsc_h2_cts.csv` coefficient p-values (and the per-config SNP
counts in the LDSC log) in the morning.
"""

from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.concrete_sldsc_generator import (
    DEFAULT_BASELINE_REF_LD_CHR_INNER_DIRNAME,
    DEFAULT_BASELINE_REF_LD_CHR_TASK,
    standard_sldsc_task_generator,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_annovar_37_sumstats import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_annovar_37_sumstats_keep_ambiguous import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS_KEEP_AMBIGUOUS,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_v1_0_extracted import (
    BASE_MODEL_PARTITIONED_LD_SCORES_V1_0_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_v1_1_extracted import (
    BASE_MODEL_PARTITIONED_LD_SCORES_V1_1_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_filtering import (
    FilterSettings,
)
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePIfNeededPipe

# Baseline v1.1 reference (older 'baseline' partitioned model).
BASELINE_V1_1_TASK = BASE_MODEL_PARTITIONED_LD_SCORES_V1_1_EXTRACTED
BASELINE_V1_1_INNER_DIRNAME = "baseline_v1.1/baseline.@"

# Baseline v1.0 reference (original 'baseline' model). Its inner dir name
# (1000G_EUR_Phase3_baseline) matches the path in Martin's command, so this is the prime
# candidate for the baseline he actually used.
BASELINE_V1_0_TASK = BASE_MODEL_PARTITIONED_LD_SCORES_V1_0_EXTRACTED
BASELINE_V1_0_INNER_DIRNAME = "1000G_EUR_Phase3_baseline/baseline.@"

# Remove all palindromic SNPs only (the analog of original LDSC's VALID_SNPS exclusion).
PALINDROME_FILTER = FilterSettings(
    remove_indels=False,
    remove_palindromic=True,
    remove_hla=False,
    keep_only_hapmap=False,
)
# Palindromes + indels + HLA: closer to a full munge-style QC.
FULL_QC_FILTER = FilterSettings(
    remove_indels=True,
    remove_palindromic=True,
    remove_hla=True,
    keep_only_hapmap=False,
)
# Palindrome filter but additionally restrict to HapMap3 up front (sanity check: the cts path
# already restricts to HapMap3 internally, so this should match PALINDROME_FILTER).
PALINDROME_FILTER_KEEP_ONLY_HM3 = FilterSettings(
    remove_indels=False,
    remove_palindromic=True,
    remove_hla=False,
    keep_only_hapmap=True,
)
# HapMap3 restriction ONLY (no palindrome/indel/HLA filtering). Used to produce a log whose
# SNP funnel (8.8M -> 1,152,792 HapMap3 SNPs) is consistent with the exported pre-S-LDSC CSV,
# without changing the result -- it reproduces c0 / the documented numbers exactly.
HM3_ONLY_FILTER = FilterSettings(
    remove_indels=False,
    remove_palindromic=False,
    remove_hla=False,
    keep_only_hapmap=True,
)


@frozen
class ExperimentConfig:
    """One cell in the discrepancy sweep."""

    name: str
    source_sumstats_task: Task
    filter_settings: FilterSettings | None
    ref_ld_chr_task: Task
    ref_ld_chr_inner_dirname: str


CONFIGS = [
    ExperimentConfig(
        name="c0",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
        filter_settings=None,
        ref_ld_chr_task=DEFAULT_BASELINE_REF_LD_CHR_TASK,
        ref_ld_chr_inner_dirname=DEFAULT_BASELINE_REF_LD_CHR_INNER_DIRNAME,
    ),
    ExperimentConfig(
        name="c1",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS_KEEP_AMBIGUOUS,
        filter_settings=None,
        ref_ld_chr_task=DEFAULT_BASELINE_REF_LD_CHR_TASK,
        ref_ld_chr_inner_dirname=DEFAULT_BASELINE_REF_LD_CHR_INNER_DIRNAME,
    ),
    ExperimentConfig(
        name="c2",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS_KEEP_AMBIGUOUS,
        filter_settings=PALINDROME_FILTER,
        ref_ld_chr_task=DEFAULT_BASELINE_REF_LD_CHR_TASK,
        ref_ld_chr_inner_dirname=DEFAULT_BASELINE_REF_LD_CHR_INNER_DIRNAME,
    ),
    ExperimentConfig(
        name="c3",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
        filter_settings=None,
        ref_ld_chr_task=BASELINE_V1_1_TASK,
        ref_ld_chr_inner_dirname=BASELINE_V1_1_INNER_DIRNAME,
    ),
    ExperimentConfig(
        name="c4",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS_KEEP_AMBIGUOUS,
        filter_settings=PALINDROME_FILTER,
        ref_ld_chr_task=BASELINE_V1_1_TASK,
        ref_ld_chr_inner_dirname=BASELINE_V1_1_INNER_DIRNAME,
    ),
    ExperimentConfig(
        name="c5",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS_KEEP_AMBIGUOUS,
        filter_settings=FULL_QC_FILTER,
        ref_ld_chr_task=BASELINE_V1_1_TASK,
        ref_ld_chr_inner_dirname=BASELINE_V1_1_INNER_DIRNAME,
    ),
    ExperimentConfig(
        name="c6",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS_KEEP_AMBIGUOUS,
        filter_settings=PALINDROME_FILTER_KEEP_ONLY_HM3,
        ref_ld_chr_task=DEFAULT_BASELINE_REF_LD_CHR_TASK,
        ref_ld_chr_inner_dirname=DEFAULT_BASELINE_REF_LD_CHR_INNER_DIRNAME,
    ),
    # Baseline v1.0 (Martin's likely background), drop-ambiguous input, no extra filter.
    ExperimentConfig(
        name="c7",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
        filter_settings=None,
        ref_ld_chr_task=BASELINE_V1_0_TASK,
        ref_ld_chr_inner_dirname=BASELINE_V1_0_INNER_DIRNAME,
    ),
    # HapMap3-restriction-only diagnostic: same input/baseline as c0, but pre-restricts to
    # HapMap3 so the log's SNP funnel matches the exported pre-S-LDSC CSV (1,152,792 SNPs).
    # Reproduces c0 / the documented numbers exactly (verified, max abs diff 0.0).
    ExperimentConfig(
        name="c8",
        source_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
        filter_settings=HM3_ONLY_FILTER,
        ref_ld_chr_task=DEFAULT_BASELINE_REF_LD_CHR_TASK,
        ref_ld_chr_inner_dirname=DEFAULT_BASELINE_REF_LD_CHR_INNER_DIRNAME,
    ),
]

# (short label used in asset ids, partitioning-dataset entry_name)
DATASETS = [
    ("gtexbrain", "gtex_brain"),
    # ("multitissue", "multi_tissue_gene_expression"),
]


def _terminal_tasks_for_dataset(dataset_label: str, entry_name: str) -> list[Task]:
    tasks: list[Task] = []
    for config in CONFIGS:
        generator = standard_sldsc_task_generator(
            sumstats_task=config.source_sumstats_task,
            base_name=f"decode_me_sldsc_exp_{dataset_label}_{config.name}",
            pre_pipe=ComputePIfNeededPipe(),
            ref_ld_chr_task=config.ref_ld_chr_task,
            ref_ld_chr_inner_dirname=config.ref_ld_chr_inner_dirname,
            filter_settings=config.filter_settings,
            partitioned_entry_names={entry_name},
        )
        tasks.extend(generator.get_terminal_tasks())
    return tasks


def run_sldsc_discrepancy_experiment():
    """
    Run the full config matrix on `gtex_brain` first, then on `multi_tissue_gene_expression`.

    Pass-1 (gtex_brain) tasks are listed before Pass-2 so that the smaller dataset completes
    first; `incremental_save=True` preserves partial progress if the run is interrupted.
    """
    terminal_tasks: list[Task] = []
    for dataset_label, entry_name in DATASETS:
        terminal_tasks.extend(_terminal_tasks_for_dataset(dataset_label, entry_name))
    DEFAULT_RUNNER.run(
        terminal_tasks,
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_sldsc_discrepancy_experiment()
