import numpy as np
from pathlib import Path

import narwhals as nw
import matplotlib.pyplot as plt
import polars as pl
import pandas as pd
import scipy.sparse
import seaborn as sns

from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import align_gwas_and_ld
from mecfs_bio.constants.gwaslab_constants import GWASLAB_NON_EFFECT_ALLELE_COL, GWASLAB_EFFECT_ALLELE_COL, \
    GWASLAB_POS_COL, GWASLAB_CHROM_COL, GWASLAB_RSID_COL
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


def go():
    # gwas_table = pl.read_parquet("assets/base_asset_store/gwas/ME_CFS/DecodeME/processed/decode_me_1_harmonize_with_ld_chr1_173000001_via_alleles.parquet")
    # ld_labels = pl.read_csv("assets/base_asset_store/reference_data/ukbb_reference_ld/chr1_173000001_17600000/raw/chr1_173000001_176000001.gz",
    #                             separator="\t")
    # pipe = CompositePipe(
    #     [
    #         RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
    #         RenameColPipe(old_name="chromosome", new_name=GWASLAB_CHROM_COL),
    #         RenameColPipe(old_name="position", new_name=GWASLAB_POS_COL),
    #         RenameColPipe(
    #             old_name="allele1",
    #             new_name=GWASLAB_NON_EFFECT_ALLELE_COL,  # See: https://github.com/omerwe/polyfun/issues/208#issuecomment-2563832487
    #         ),
    #         RenameColPipe(old_name="allele2", new_name=GWASLAB_EFFECT_ALLELE_COL),
    #     ],
    # )
    # ld_labels=pipe.process(nw.from_native(ld_labels).lazy()).collect().to_polars()
    #
    # ld_matrix = scipy.sparse.load_npz("assets/base_asset_store/reference_data/ukbb_reference_ld/chr1_173000001_17600000/raw/chr1_173000001_176000001.npz")
    #
    fig, axes= plt.subplots(1, 1,figsize=(12,12))
    # gwas, ld, ld_mat=align_gwas_and_ld(
    #     gwas=gwas_table,
    #     ld_labels=ld_labels,
    #     ld_matrix_sparse=ld_matrix,
    # )
    # gwas.write_parquet("data/gwas.parquet")
    # ld.write_parquet("data/ld.parquet")
    # np.savez("data/corr_mat.npz", ld_mat)
    # gwas = pl.read_parquet("data/gwas.parquet")
    # ld = pl.read_parquet("data/ld.parquet")
    loaded=  np.load("data/corr_mat.npz")
    ld_mat =loaded["arr_0"]
    ld_mat = ld_mat+ld_mat.T
    ld_mat_squared = abs(ld_mat)
    ax = sns.heatmap(
        ld_mat_squared,
        # annot=True,  # Add correlation values as annotations
        # fmt='.2f',  # Format annotations to two decimal places
        cmap='coolwarm',  # Use a diverging colormap
        # vmin=-1,  # Ensure color scale covers full range from -1 to 1
        vmax=1,
        center=0,  # Center the colormap at 0
        square=True,  # Make cells square shaped
        # linewidths=.5,  # Add lines to divide cells
        cbar_kws={"shrink": .75} , # Adjust color bar size
        ax=axes,
        xticklabels=False, yticklabels=False

    )
    ax.tick_params(left=False, bottom=False)
    fig.show()
    import pdb; pdb.set_trace()
    write_plots_to_dir(
       Path( "data"),
        {
            "my_corr_plot_2":fig
        }

    )
    import pdb; pdb.set_trace()


if __name__ == '__main__':
    go()