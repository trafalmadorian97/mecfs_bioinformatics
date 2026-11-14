"""
For each display-table hit of the exploratory DecodeME x PPP coloc run, checks whether the
trait/protein effect-sign agreement holds across the shared signal, not only at the protein lead:
among aligned variants with protein p < 1e-6 and trait p < 1e-4, the fraction whose trait and
protein z-scores (same effect allele) agree in sign, plus the sign agreement at the trait lead.

Usage: pixi r python experiments/claude/protein_coloc/sign_consistency_probe.py
"""

from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_coloc import (
    DECODE_ME_PPP_COLOC_EXPLORATORY_2,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    two_sided_mlog10p,
)
from mecfs_bio.constants.protein_coloc_constants import (
    BETA_PROTEIN_COL,
    BETA_TRAIT_COL,
    DISPLAY_PROTEIN_COL,
    DISPLAY_LOCUS_COL,
    DISPLAY_ASSAY_ID_COL,
    DISPLAY_SAME_SIGN_COL,
    SE_PROTEIN_COL,
    SE_TRAIT_COL,
    WINDOWS_DIRNAME,
    pair_file_stem,
)

if __name__ == "__main__":
    group = DECODE_ME_PPP_COLOC_EXPLORATORY_2
    coloc_dir = Path(DEFAULT_RUNNER.meta_to_path(group.coloc_task.meta))
    display = pl.read_parquet(DEFAULT_RUNNER.meta_to_path(group.display_task.meta))
    rows = []
    for gene, protein, locus, same_sign in display.select(
        DISPLAY_PROTEIN_COL, DISPLAY_ASSAY_ID_COL, DISPLAY_LOCUS_COL, DISPLAY_SAME_SIGN_COL
    ).iter_rows():
        w = pl.read_parquet(
            coloc_dir / WINDOWS_DIRNAME / f"{pair_file_stem(locus, protein)}.parquet"
        )
        zt = (w[BETA_TRAIT_COL] / w[SE_TRAIT_COL]).to_numpy()
        zp = (w[BETA_PROTEIN_COL] / w[SE_PROTEIN_COL]).to_numpy()
        shared = (two_sided_mlog10p(zp) > 6) & (two_sided_mlog10p(zt) > 4)
        trait_lead = int(np.argmax(np.abs(zt)))
        rows.append(
            {
                "gene": gene,
                "same_sign_at_protein_lead": same_sign,
                "same_sign_at_trait_lead": bool(np.sign(zt[trait_lead]) == np.sign(zp[trait_lead])),
                "n_shared": int(shared.sum()),
                "frac_shared_agree": float((np.sign(zt[shared]) == np.sign(zp[shared])).mean())
                if shared.any()
                else None,
            }
        )
    with pl.Config(tbl_cols=10, tbl_width_chars=200):
        print(pl.DataFrame(rows))
