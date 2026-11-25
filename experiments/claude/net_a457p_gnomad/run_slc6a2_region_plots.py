"""Build the SLC6A2 region plots for DecodeME GWAS-1 and Keaton DBP.

Run from the main checkout (so the asset store and its config are found):
pixi r python experiments/claude/net_a457p_gnomad/run_slc6a2_region_plots.py
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_region_plot_slc6a2_locus_37 import (
    KEATON_DBP_SLC6A2_REGION_PLOT_37,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_region_plot_slc6a2_locus_37 import (
    DECODE_ME_SLC6A2_REGION_PLOT_37,
)

if __name__ == "__main__":
    DEFAULT_RUNNER.run(
        [DECODE_ME_SLC6A2_REGION_PLOT_37, KEATON_DBP_SLC6A2_REGION_PLOT_37],
        incremental_save=True,
    )
