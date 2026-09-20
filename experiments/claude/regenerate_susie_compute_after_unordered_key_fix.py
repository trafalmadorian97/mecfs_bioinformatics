"""Force-rebuild the SUSIE compute tasks that feed the SUSIE figures, so their
outputs reflect the retired-unordered-allele-key correction (PR #1181).

Why this script exists
----------------------
The verifying-trace rebuilder hashes asset files, not Python task code, so the
allele-key correction -- which changed the SUSIE gwas harmonization pipe, the
exact (chrom, pos, ea, nea) prior join, the exact-key annotation dedup, and the
annotation ridge-weight join -- does not by itself invalidate the cached SUSIE
assets. This script force-rebuilds the changed compute tasks (and, via
must_rebuild_transitive, every task downstream of them within the target graph)
so the corrected outputs land in the asset store.

It deliberately rebuilds only the compute tasks (SUSIE fits, the polyfun
explain contrasts, and the two annotation reference assets), not the figures.
The figures are re-plotted and copied into the figure directory afterwards by
regenerate_susie_figures_after_unordered_key_fix.py, which reads these
now-corrected cached outputs.

Run:
    pixi r python \
      experiments/claude/regenerate_susie_compute_after_unordered_key_fix.py \
      2>&1 | tee experiments/claude/logs/regenerate_susie_compute_after_unordered_key_fix.log

Dry run (print the target and force-rebuild asset ids, build nothing):
    pixi r python \
      experiments/claude/regenerate_susie_compute_after_unordered_key_fix.py --dry-run
"""

import sys

import structlog

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr1_174_128_548 import (
    POLYFUN_EXPLAIN_CHR1_174,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr15_54_925_638 import (
    POLYFUN_EXPLAIN_CHR15_54,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr17_50_237_377 import (
    POLYFUN_EXPLAIN_CHR17_50,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_explainability.susie_explain_decode_me_37_chr20_47_653_230 import (
    POLYFUN_EXPLAIN_CHR20_47,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr1_174_128_548_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_97_505_620_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr15_54_925_638_locus_plalindindromes import (
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr17_50_237_377_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr20_47_653_230_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.annotation_ridge_weights import (
    BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.asset_generator.fine_mapping_asset_generator import (
    BroadFineMapTaskGroup,
)
from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    PolyfunExplainOuterGroup,
)
from mecfs_bio.build_system.task.base_task import Task

logger = structlog.get_logger()

# The five with-palindromes fine-mapping loci. Each exposes four SUSIE fits
# (L=10, L=10 strict, L=1, L=2). These runs carry no polyfun prior, so their
# only changed input is the harmonization gwas pipe now applied inside the SUSIE
# task.
WITH_PALINDROME_GROUPS: list[BroadFineMapTaskGroup] = [
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES,
]

# The four polyfun-explainability loci. Each outer group holds four matched
# uniform/polyfun SUSIE pairs plus a contrast per pair; the contrasts also
# depend on the two changed annotation reference assets.
POLYFUN_EXPLAIN_GROUPS: list[PolyfunExplainOuterGroup] = [
    POLYFUN_EXPLAIN_CHR1_174,
    POLYFUN_EXPLAIN_CHR15_54,
    POLYFUN_EXPLAIN_CHR17_50,
    POLYFUN_EXPLAIN_CHR20_47,
]


def _with_palindrome_susie_tasks() -> list[Task]:
    tasks: list[Task] = []
    for group in WITH_PALINDROME_GROUPS:
        tasks += [
            group.susie_finemap_task,
            group.susie_finemap_strict_task,
            group.susie_finemap_1_credible_set_task,
            group.susie_finemap_2_credible_set_task,
        ]
    return tasks


def _polyfun_compute_tasks() -> list[Task]:
    tasks: list[Task] = []
    for outer in POLYFUN_EXPLAIN_GROUPS:
        for group in outer.groups:
            tasks += [group.susie_uniform, group.susie_polyfun, group.contrast]
    return tasks


def _dedup_by_asset_id(tasks: list[Task]) -> list[Task]:
    """Two loci can share a compute task (e.g. the annotation assets appear once
    per contrast). Keep the first occurrence of each asset id."""
    seen: set[str] = set()
    unique: list[Task] = []
    for task in tasks:
        key = str(task.asset_id)
        if key not in seen:
            seen.add(key)
            unique.append(task)
    return unique


# The two annotation reference assets whose builders changed (exact-key dedup and
# exact-key ridge join). They feed the polyfun explain contrasts, not the SUSIE
# fits themselves.
CHANGED_ANNOTATION_ASSETS: list[Task] = [
    BASELINE_LF_ANNOTATION_MATRIX,
    BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS,
]

MUST_REBUILD: list[Task] = _dedup_by_asset_id(
    CHANGED_ANNOTATION_ASSETS
    + _with_palindrome_susie_tasks()
    + _polyfun_compute_tasks()
)

# Targeting the compute tasks (not the figures) keeps the plot/table tasks out of
# this run; they are re-plotted afterwards by the figure script. Every changed
# compute task is both a target and a forced rebuild.
TARGETS: list[Task] = MUST_REBUILD


def _print_plan() -> None:
    logger.info(f"{len(TARGETS)} target compute tasks; forcing rebuild of all")
    for task in TARGETS:
        logger.info(f"  rebuild: {task.asset_id}")


def regenerate_susie_compute() -> None:
    _print_plan()
    DEFAULT_RUNNER.run(
        targets=TARGETS,
        incremental_save=True,
        must_rebuild_transitive=MUST_REBUILD,
    )


if __name__ == "__main__":
    if "--dry-run" in sys.argv:
        _print_plan()
    else:
        regenerate_susie_compute()
