"""
Force-rebuilds the protein info table of the exploratory DecodeME x UKB-PPP coloc run (code changes
do not invalidate the build cache, and the table's columns changed). Only the info task is rebuilt;
the coloc and display tasks come from the cache.

Usage: pixi r python experiments/claude/protein_coloc/rebuild_decode_me_exploratory_2_protein_info.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO,
)

if __name__ == "__main__":
    DEFAULT_RUNNER.run(
        [DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO],
        incremental_save=True,
        must_rebuild_transitive=[DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO],
    )
    pl.Config.set_tbl_cols(20)
    pl.Config.set_tbl_width_chars(250)
    pl.Config.set_fmt_str_lengths(60)
    path = DEFAULT_RUNNER.meta_to_path(DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO.meta)
    print(path)
    print(pl.read_parquet(path))
