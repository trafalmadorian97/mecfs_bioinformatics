"""Re-plot and re-export the SUSIE-related figures after the retired-unordered-
allele-key correction (PR #1181).

This is the second half of the refresh. First run
regenerate_susie_compute_after_unordered_key_fix.py, which force-rebuilds the
SUSIE fits, polyfun explain contrasts, and the two annotation reference assets
so their corrected outputs are cached. Then run this script, which calls
regenerate_figures on exactly the SUSIE-related figures in ALL_FIGURE_TASKS.

regenerate_figures force-rebuilds the plotting/table tasks (reading the now
corrected cached compute outputs) and copies the results into the figure
directory. It does not push blobs to the figures release; do that separately if
the refreshed figures should ship.

Run:
    pixi r python \
      experiments/claude/regenerate_susie_figures_after_unordered_key_fix.py \
      2>&1 | tee experiments/claude/logs/regenerate_susie_figures_after_unordered_key_fix.log

Dry run (print the SUSIE figure asset ids, build nothing):
    pixi r python \
      experiments/claude/regenerate_susie_figures_after_unordered_key_fix.py --dry-run
"""

import sys

import structlog

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
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures

logger = structlog.get_logger()

# The SUSIE-related figure tasks, grouped exactly as they appear in
# ALL_FIGURE_TASKS: the with-palindromes fine-mapping figures + credible-set
# markdown tables, and the polyfun-explainability figures + tables. Built from
# the generator objects (not re-derived by hand) and filtered against
# ALL_FIGURE_TASKS below, so a task that is not actually a published figure is
# dropped rather than silently re-exported.
_WITH_PALINDROME_FIGURES: list[Task] = [
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.susie_stackplot_task,
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.upset_plot_task,
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.susie_base_credible_set_markdown_table,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.susie_stackplot_task,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.susie_finemap_2_credible_set_plot,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.upset_plot_task,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.susie_base_credible_set_markdown_table,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.susie_2_credible_set_markdown_table,
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES.susie_stackplot_task,
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES.upset_plot_task,
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES.susie_base_credible_set_markdown_table,
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.susie_stackplot_task,
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.upset_plot_task,
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.upset_plot_task_pip001,
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.susie_base_credible_set_markdown_table,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_stackplot_task,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_finemap_2_credible_set_plot,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_finemap_strict_plot,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.upset_plot_task,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_base_credible_set_markdown_table,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_strict_credible_set_markdown_table,
]


def _polyfun_figures() -> list[Task]:
    figures: list[Task] = [
        POLYFUN_EXPLAIN_CHR1_174.upset_all_polyfun,
        POLYFUN_EXPLAIN_CHR1_174.upset_cs50_polyfun,
        POLYFUN_EXPLAIN_CHR15_54.upset_all_polyfun,
        POLYFUN_EXPLAIN_CHR15_54.upset_cs50_polyfun,
        POLYFUN_EXPLAIN_CHR17_50.upset_all_polyfun,
        POLYFUN_EXPLAIN_CHR17_50.upset_cs50_polyfun,
        POLYFUN_EXPLAIN_CHR20_47.upset_all_polyfun,
        POLYFUN_EXPLAIN_CHR20_47.upset_cs50_polyfun,
    ]
    # chr1/15/17 publish only the l10 pair's plot + two tables; chr20 also
    # publishes the l1 pair's.
    for outer in (
        POLYFUN_EXPLAIN_CHR1_174,
        POLYFUN_EXPLAIN_CHR15_54,
        POLYFUN_EXPLAIN_CHR17_50,
        POLYFUN_EXPLAIN_CHR20_47,
    ):
        group = outer.groups_by_label["l10"]
        figures += [
            group.plot_svg,
            group.detailed_table,
            group.per_variant_annotation_table,
        ]
    l1_chr20 = POLYFUN_EXPLAIN_CHR20_47.groups_by_label["l1"]
    figures += [
        l1_chr20.plot_svg,
        l1_chr20.detailed_table,
        l1_chr20.per_variant_annotation_table,
    ]
    return figures


# Keep only the ones actually published in ALL_FIGURE_TASKS.
SUSIE_FIGURE_TASKS: list[Task] = [
    task
    for task in _WITH_PALINDROME_FIGURES + _polyfun_figures()
    if task in ALL_FIGURE_TASKS
]


def _print_plan() -> None:
    logger.info(f"{len(SUSIE_FIGURE_TASKS)} SUSIE figure tasks to regenerate")
    for task in SUSIE_FIGURE_TASKS:
        logger.info(f"  figure: {task.asset_id}")


if __name__ == "__main__":
    if "--dry-run" in sys.argv:
        _print_plan()
    else:
        _print_plan()
        regenerate_figures(SUSIE_FIGURE_TASKS)
