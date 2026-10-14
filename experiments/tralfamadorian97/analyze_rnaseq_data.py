import matplotlib.pyplot as plt
from pathlib import Path

import numpy as np
import seaborn as sns
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
#
#
# def analyze_matrix():
#     pca = PCA(n_components=10)
import pandas as pd

from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


def get_prep_for_magma_frame()-> pd.DataFrame:
    return pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/processed/gtex_v10_rna_seq_median_tissue_expression_prep_for_magma.tsv",
                sep="\t",)

def scale_frame(df: pd.DataFrame) -> pd.DataFrame:
    scaler= StandardScaler()
    scaled = scaler.fit_transform(df)
    return pd.DataFrame(scaled, columns=df.columns)

def drop_zero_genes(df: pd.DataFrame) -> pd.DataFrame:
    zero_genes= (df==0).all(axis=1)
    return df[~zero_genes]


def get_pca_components_frame(rna_seq_df: pd.DataFrame) -> pd.DataFrame:
    pca = PCA(n_components=2)
    pca.fit(rna_seq_df)
    comp_df= pd.DataFrame(pca.components_, columns=pd.Index(rna_seq_df.columns,name="tissue"), index=pd.Index(["comp1", "comp2"], name="component"))
    return comp_df


def go():
    df= get_prep_for_magma_frame().drop("Name", axis=1)
    df= drop_zero_genes(df)
    # scaled_df = scale_frame(df)
    comp_df = get_pca_components_frame(df)
    comp_df_t = comp_df.transpose().reset_index()
    plots = {}
    plots["tissue_pca_components_no_scaling"]= px.scatter(comp_df_t, x="comp1", y="comp2", text="tissue")
    write_plots_to_dir(Path("output/rna_figs"), plots)
    import pdb; pdb.set_trace()
    print("yo")


def plot_median_hist():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                       sep="\t",header=2)
    df= drop_zero_genes(df)
    df=  df.drop(["Name","Description"], axis=1)
    flat_df = pd.DataFrame(df.values.reshape(-1), columns = ["val"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    axs.hist(flat_df.values, bins=100)
    # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["hist_of_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)

def plot_median_ecdf():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df= drop_zero_genes(df)
    import pdb; pdb.set_trace()
    df=  df.drop(["Name","Description"], axis=1)
    flat_df = pd.DataFrame(df.values.reshape(-1), columns = ["val"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    sns.ecdfplot(flat_df, x="val", ax=axs)
    # axs.hist(flat_df.values, bins=100)
    # # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["ecdf_of_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)


def plot_median_ecdf_log_nonzero():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df= drop_zero_genes(df)
    df=  df.drop(["Name","Description"], axis=1)
    flat = df.values.reshape(-1)
    nonzero_flat = flat[flat>0]
    log_nonzero_flat = np.log(nonzero_flat)
    import pdb; pdb.set_trace()
    flat_df = pd.DataFrame(log_nonzero_flat, columns = ["log_tpm"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    sns.ecdfplot(flat_df, x="log_tpm", ax=axs)
    # axs.hist(flat_df.values, bins=100)
    # # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["ecdf_of_log_of_nonzero_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)


def plot_median_hist_log_nonzero():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df= drop_zero_genes(df)
    df=  df.drop(["Name","Description"], axis=1)
    flat = df.values.reshape(-1)
    nonzero_flat = flat[flat>0]
    log_nonzero_flat = np.log(nonzero_flat)
    import pdb; pdb.set_trace()
    flat_df = pd.DataFrame(log_nonzero_flat, columns = ["log_tpm"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    sns.histplot(flat_df, x="log_tpm", ax=axs)
    # axs.hist(flat_df.values, bins=100)
    # # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["hist_of_log_of_nonzero_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)

def plot_median_ecdf_shifted_log_all():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df= drop_zero_genes(df)
    df=  df.drop(["Name","Description"], axis=1)
    flat = df.values.reshape(-1)
    shifted_log_flat = np.log(flat+1)
    import pdb; pdb.set_trace()
    flat_df = pd.DataFrame(shifted_log_flat, columns = ["shifted_log_tpm"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    sns.ecdfplot(flat_df, x="shifted_log_tpm", ax=axs)
    # axs.hist(flat_df.values, bins=100)
    # # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["ecdf_of_shifted_log_of_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)


def plot_median_ecdf_shifted_log_v2_all():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df= drop_zero_genes(df)
    df=  df.drop(["Name","Description"], axis=1)
    flat = df.values.reshape(-1)
    basic_flat_df = pd.DataFrame(flat, columns = ["tpm"])
    nonzero_flat_df = pd.DataFrame(flat[flat>0], columns = ["tpm"])
    import pdb; pdb.set_trace()
    min_pos = np.min(flat[flat>0])
    shifted_log_flat = np.log(flat+min_pos)
    flat_df = pd.DataFrame(shifted_log_flat, columns = ["shifted_log_tpm_v2"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    sns.ecdfplot(flat_df, x="shifted_log_tpm_v2", ax=axs)
    # axs.hist(flat_df.values, bins=100)
    # # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["ecdf_of_shifted_log_v2_of_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)


def plot_median_ecdf_shifted_log_v3_all():
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df= drop_zero_genes(df)
    df=  df.drop(["Name","Description"], axis=1)
    flat = df.values.reshape(-1)
    basic_flat_df = pd.DataFrame(flat, columns = ["tpm"])
    nonzero_flat_df = pd.DataFrame(flat[flat>0], columns = ["tpm"])
    psuedocount = 0.01
    shifted_log_flat = np.log(flat+psuedocount)
    flat_df = pd.DataFrame(shifted_log_flat, columns = ["shifted_log_tpm_v3"])
    fig, axs = plt.subplots(nrows=1, ncols=1)
    sns.ecdfplot(flat_df, x="shifted_log_tpm_v3", ax=axs)
    # axs.hist(flat_df.values, bins=100)
    # # ecdf = px.ecdf(flat_df, x="val")
    plots = {}
    plots["ecdf_of_shifted_log_v3_of_median_values"] = fig
    write_plots_to_dir(Path("output/rna_figs"), plots)



def try_pca_with_adjusted_psuedocount():
    # df= get_prep_for_magma_frame().drop("Name", axis=1)
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                     sep="\t",header=2)
    df=  df.drop(["Name","Description"], axis=1)
    df= drop_zero_genes(df)
    df_pseudo = np.log(df+0.01)
    import pdb; pdb.set_trace()
    # scaled_df = scale_frame(df)
    comp_df = get_pca_components_frame(df_pseudo)
    comp_df_t = comp_df.transpose().reset_index()
    plots = {}
    plots["tissue_pca_components_001_psuedocount"]= px.scatter(comp_df_t, x="comp1", y="comp2", text="tissue")
    write_plots_to_dir(Path("output/rna_figs"), plots)
    import pdb; pdb.set_trace()
    print("yo")



def try_pca_with_1_filter():
    # df= get_prep_for_magma_frame().drop("Name", axis=1)
    df = pd.read_csv("assets/base_asset_store/reference_data/rna_seq_data/gtex/extracted/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",sep="\t",header=2)
    df=  df.drop(["Name","Description"], axis=1)
    df= drop_zero_genes(df)
    df = df.loc[(df>=1).any(axis=1),: ] # drop genes that do not exceed 1 TPM in some tissue
    df_pseudo = np.log(df+1)
    # scaled_df = scale_frame(df)
    comp_df = get_pca_components_frame(df_pseudo)
    comp_df_t = comp_df.transpose().reset_index()
    plot  =px.scatter(comp_df_t, x="comp1", y="comp2", text="tissue")
    sel = comp_df_t.loc[:,["comp1", "comp2","tissue"]]
    print(sel.to_json( index=False))
    import pdb; pdb.set_trace()
    plots = {}
    plots["tissue_pca_components_with_1_filter"]=  plot
    json_str = plot.to_json()
    print(json_str)
    write_plots_to_dir(Path("output/rna_figs"), plots)
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == "__main__":
    # try_pca_with_adjusted_psuedocount()
    try_pca_with_1_filter()
    # plot_median_ecdf_shifted_log_v3_all()
    # plot_median_ecdf_log_nonzero():W

    # plot_median_hist_log_nonzero()
    # plot_median_ecdf_shifted_log_all()
    # plot_median_hist()
    # plot_median_ecdf()
    # go()
