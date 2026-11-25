"""Full coloc.abf posteriors (H0-H4) of the AD-run pairs for the Western AD aptamers that
tested below PP.H4 0.8 (PRSS8, LILRB1, CA12)."""

import polars as pl

from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_alzheimers_csf_coloc import (
    BELLENGUEZ_ALZHEIMERS_CSF_COLOC,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.constants.protein_coloc_constants import (
    RAW_RESULTS_FILENAME,
)

IDS = ["15513-108", "5090-49", "19347-37"]

out_dir = DEFAULT_RUNNER.meta_to_path(BELLENGUEZ_ALZHEIMERS_CSF_COLOC.coloc_task.meta)
results = pl.read_parquet(out_dir / RAW_RESULTS_FILENAME)
print(results.columns)
pl.Config.set_tbl_cols(30)
pl.Config.set_tbl_width_chars(250)


cols = [c for c in results.columns if c.startswith("PP_H") or "NSNP" in c or "N_SNP" in c]
print(results.filter(pl.col("PROTEIN_ID").is_in(IDS)).select("GENE", *cols))
