from pathlib import Path

import pandas as pd

from mecfs_bio.constants.gwaslab_constants import GWASLAB_BETA_COL, GWASLAB_SE_COL

file_sent_to_martin_path = Path("martin_exchange_5/hapmap_filtered_sumstats.csv")
martin_munged_path = Path("gwas_1.regenie.filtered.rsids.munged.txt.sumstats.gz")


def _prep_munged_df_martin() -> pd.DataFrame:
    df_martin = pd.read_csv(martin_munged_path, sep="\s+")
    df_martin = df_martin.dropna(subset=["Z"])
    return df_martin

def _prep_df_peter() ->pd.DataFrame:

    df_peter = pd.read_csv(file_sent_to_martin_path)
    return df_peter


def intersect_rsids():
    """
    Output:
    Number of SNPs (Peter): 1152792
    Number of SNPs (Martin): 1157905
    Number of common SNPs 1152587
    Number of Peter-Unique: SNPs 205
    Number of martin-Unique: SNPs 5318
    """
    df_peter = _prep_df_peter()
    peter_rsids=set(df_peter["rsID"].tolist())

    df_martin = _prep_munged_df_martin()
    martin_rsids = set(df_martin["SNP"].tolist())

    print(f"Number of SNPs (Peter): {len(peter_rsids)}")
    print(f"Number of SNPs (Martin): {len(martin_rsids)}")

    print(f"Number of common SNPs: {len(peter_rsids & martin_rsids)}")

    print(f"Number of Peter-Unique SNPs: {len(peter_rsids - martin_rsids)}")

    print(f"Number of Martin-Unique SNPs: {len(martin_rsids-peter_rsids)}")

def add_zscore_to_df_peter(df_peter:pd.DataFrame,suffix:str)-> pd.DataFrame:
    df_peter = df_peter.copy()
    df_peter["Z"+suffix] = df_peter[GWASLAB_BETA_COL]/df_peter[GWASLAB_SE_COL]
    return df_peter


def check_zscore_agreement():
    df_peter = _prep_df_peter()
    df_peter= add_zscore_to_df_peter(df_peter, suffix="_peter")
    df_martin = _prep_munged_df_martin()
    joined = df_martin.merge(df_peter, left_on="SNP", right_on="rsID")
    joined["z_diff"] = abs( joined["Z"] - joined["Z_peter"]   )
    print("Absolute z diff stats:")
    print(joined["z_diff"].describe(percentiles=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]))
    print("Relative z diff stats:")
    joined = joined.loc[abs(joined["Z"])>0]
    joined["rel_z_diff"] = abs( (joined["Z"] - joined["Z_peter"] )/abs(joined[["Z","Z_peter"]].max(axis=1))  )
    print(joined["rel_z_diff"].describe(percentiles=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]))
    import pdb; pdb.set_trace()
    print("yo")






if __name__ == "__main__":
    check_zscore_agreement()
    # intersect_rsids()

