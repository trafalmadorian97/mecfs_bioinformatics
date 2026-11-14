"""
Force-rebuilds the docs-facing PPP coloc display tables (code changes do not invalidate the build
cache, the display columns were renamed: gene -> protein, protein -> assay_id, and the DecodeME
stack plots label the trait row "DecodeME GWAS-1"), then copies the
display tables, the DecodeME stack plots and the DecodeME protein info table to the figure
directory. The stack plots and protein info table rebuild because their display-table dependency
changed; the coloc results come from the cache.

Usage: pixi r python experiments/claude/protein_coloc/rebuild_coloc_display_assets.py 2>&1 | tee experiments/claude/protein_coloc/rebuild_coloc_display_assets.log
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_ppp_coloc import (
    WILLER_LDL_PPP_COLOC,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC_EXPLORATORY_2,
    DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
    DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS,
)
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures

if __name__ == "__main__":
    display_tasks = [
        DECODE_ME_PPP_COLOC_EXPLORATORY_2.display_task,
        WILLER_LDL_PPP_COLOC.display_task,
    ]
    figure_tasks = [
        *display_tasks,
        DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS.png_task,
        DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
    ]
    DEFAULT_RUNNER.run(
        figure_tasks, incremental_save=True, must_rebuild_transitive=display_tasks
    )
    regenerate_figures(figure_tasks)
