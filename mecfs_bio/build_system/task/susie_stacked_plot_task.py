import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
from matplotlib import gridspec

import polars as pl
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL


def plot_locus_tracks_matplotlib(
    gwas_df: pl.DataFrame,
    susie_df: pl.DataFrame,
    ld_r2: np.ndarray,
    gene_df: pl.DataFrame,
    *,
    gwas_pos_col: str = "pos",
    gwas_p_col: str = "p",
    susie_pos_col: str = "pos",
    susie_pip_col: str = "pip",
    susie_cs_col: str = "cs",
    gene_start_col: str = "start",
    gene_end_col: str = "end",
    gene_name_col: str = "name",
    gene_strand_col: str = "strand",
    start_bp: int | None = None,
    end_bp: int | None = None,
    title: str = "Locus tracks",
    ld_vmin: float = 0.0,
    ld_vmax: float = 1.0,
):
    start_dp = gene_df[gwas_pos_col].min()
    end_bp =  gene_df[gwas_pos_col].max()
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(
        nrows=4, ncols=1,
        height_ratios=[2.2, 1.2, 2.2, 1.3],
        hspace=0.08
    )
