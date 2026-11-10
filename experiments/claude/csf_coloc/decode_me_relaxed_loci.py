"""
Count DecodeME GWAS-1 trait loci (as the coloc task defines them: greedy leads, +-500 kb, merged,
MHC excluded) at lead p thresholds 5e-8, 1e-6 and 1e-5, to size the exploratory runs' pair
budget.

Usage: pixi r python experiments/claude/csf_coloc/decode_me_relaxed_loci.py
"""

import polars as pl
from scipy.stats import norm

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_38_harmonized import (
    DECODE_ME_GWAS_1_BUILD_38_HARMONIZED as TRAIT,
)
from mecfs_bio.build_system.task.protein_coloc.trait_loci import find_trait_loci
from mecfs_bio.constants.genomic_coordinate_constants import extended_mhc_interval
from mecfs_bio.constants.gwaslab_constants import GWASLAB_BETA_COL, GWASLAB_SE_COL
from mecfs_bio.constants.protein_coloc_constants import STATUS_COL

THRESHOLDS = (5e-8, 1e-6, 1e-5)


def main() -> None:
    trait = pl.scan_parquet(DEFAULT_RUNNER.meta_to_path(TRAIT.meta))
    for p in THRESHOLDS:
        z = float(norm.isf(p / 2)) * (1 - 1e-6)
        sig = trait.filter(
            (pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)).abs() >= z
        ).collect()
        loci = find_trait_loci(
            sig,
            p_threshold=p,
            half_window_bp=500_000,
            exclude=[extended_mhc_interval("38")],
        )
        print(
            f"p<{p:g}: significant variants {sig.height}, loci {loci.height}, "
            f"by status {dict(loci.group_by(STATUS_COL).len().rows())}"
        )


if __name__ == "__main__":
    main()
