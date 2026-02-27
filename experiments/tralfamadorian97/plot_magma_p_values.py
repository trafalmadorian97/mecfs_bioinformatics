from pathlib import Path

import numpy as np
import plotly.express as px
import pandas as pd

from mecfs_bio.constants.gwaslab_constants import GWASLAB_MLOG10P_COL
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


def plot_p_values():
    df = pd.read_csv("assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis_results/magma/decode_me_gwas_1_build_37_magma_ensemble_specific_tissue_gene_covar_analysis/gene_set_analysis_output.gsa.out",
                     header=5, sep="\s+" )

    num_tests = len(df)
    df[GWASLAB_MLOG10P_COL] = -np.log10(df["P"])
    print(f"Found {num_tests} tests")
    thresh = 0.05/num_tests
    print(f"Bonferroni thresh: {thresh}")
    trunc_df = df.sort_values(by=['P'], ascending=True).iloc[:20]
    trunc_df["Significance"] = float("nan")
    trunc_df["Significance"] = "Significant"
    trunc_df["Significance"] = trunc_df["Significance"].where(trunc_df["P"]<=thresh, "Not Significant")
    plot = px.bar(trunc_df, x="FULL_NAME", y=GWASLAB_MLOG10P_COL, color="Significance")
    plots ={"decode_me_gtex_tissue":plot}
    write_plots_to_dir(
        Path("output")/"gwas_1"/"magma", plots
    )

if __name__ == "__main__":
    plot_p_values()