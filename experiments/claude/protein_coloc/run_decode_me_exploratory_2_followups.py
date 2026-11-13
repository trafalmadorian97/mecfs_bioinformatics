"""
Builds the follow-ups of the exploratory DecodeME x UKB-PPP coloc run: the display table (forced,
since code changes do not invalidate the build cache and it gained columns), the per-locus stack
plots with their SVG and PNG copy-outs, and the protein annotations. The coloc task itself is not rerun.

Usage: pixi r python experiments/claude/protein_coloc/run_decode_me_exploratory_2_followups.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC_EXPLORATORY_2,
    DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
    DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS,
)

if __name__ == "__main__":
    DEFAULT_RUNNER.run(
        [
            *DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS.terminal_tasks(),
            DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
        ],
        incremental_save=True,
        must_rebuild_transitive=[DECODE_ME_PPP_COLOC_EXPLORATORY_2.display_task],
    )
    pl.Config.set_tbl_cols(20)
    pl.Config.set_tbl_width_chars(250)
    pl.Config.set_fmt_str_lengths(120)
    for task in (
        DECODE_ME_PPP_COLOC_EXPLORATORY_2.display_task,
        DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
    ):
        print(task.asset_id)
        print(pl.read_parquet(DEFAULT_RUNNER.meta_to_path(task.meta)))
    for plot_task in (
        DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS.plot_task,
        DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS.svg_task,
        DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS.png_task,
    ):
        plots_dir = DEFAULT_RUNNER.meta_to_path(plot_task.meta)
        print(plots_dir)
        print(sorted(p.name for p in plots_dir.iterdir()))
