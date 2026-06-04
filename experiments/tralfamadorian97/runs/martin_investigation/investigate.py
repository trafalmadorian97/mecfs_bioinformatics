from pathlib import Path

import pandas as pd

from mecfs_bio.constants.gwaslab_constants import GWASLAB_BETA_COL, GWASLAB_SE_COL, GWASLAB_EFFECT_ALLELE_COL, \
    GWASLAB_NON_EFFECT_ALLELE_COL

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

def add_abs_zscore_to_df_peter(df_peter:pd.DataFrame,suffix:str)-> pd.DataFrame:
    df_peter = df_peter.copy()
    df_peter["Z"+suffix] = df_peter[GWASLAB_BETA_COL]/df_peter[GWASLAB_SE_COL]
    df_peter["Z_abs"+suffix] = abs(df_peter[GWASLAB_BETA_COL]/df_peter[GWASLAB_SE_COL])
    return df_peter

def _add_z_abs(df: pd.DataFrame)-> pd.DataFrame:
    df = df.copy()
    df["Z_abs"]= abs(df["Z"])
    return df


def _get_value_counts():
    """
    passes
    """
    df_peter = _prep_df_peter()
    df_martin = _prep_munged_df_martin()
    vc_peter = df_peter["rsID"].value_counts()
    vc_martin  =df_martin.value_counts()
    assert vc_peter.max()==1
    assert vc_martin.max()==1
    print("passes")


def _get_allele_disagreement():
    """
    Forward matches 3025
    Reverse matches 1149561
    Instances of variants with the same rsID but different alleles 1
    """
    df_peter = _prep_df_peter()
    df_martin = _prep_munged_df_martin()
    joined = df_martin.merge(df_peter, left_on="SNP", right_on="rsID")
    print(f"Joined DF has {len(joined)} rows")
    # import pdb; pdb.set_trace()
    match_forward =joined.loc[ (joined[GWASLAB_EFFECT_ALLELE_COL]==joined["A1"]) &   (joined[GWASLAB_NON_EFFECT_ALLELE_COL]==joined["A2"])  ]
    match_reverse =joined.loc[ (joined[GWASLAB_EFFECT_ALLELE_COL]==joined["A2"]) &   (joined[GWASLAB_NON_EFFECT_ALLELE_COL]==joined["A1"])  ]
    match_total = pd.concat([match_forward, match_reverse], axis=0)
    disagreement = joined.loc[~joined["rsID"].isin(match_total["rsID"]) ]
    print(f"Forward matches {len(match_forward)}")
    print(f"Reverse matches {len(match_reverse)}")
    print(f"Instances of variants with the same rsID but different alleles {len(disagreement)}")




def check_zscore_agreement():
    """
    Absolute z diff stats:
    count    1.152587e+06
    mean     2.505455e-04
    std      5.513735e-04
    min      0.000000e+00
    10%      4.978917e-05
    20%      1.000721e-04
    30%      1.499271e-04
    40%      2.000814e-04
    50%      2.501022e-04
    60%      3.001047e-04
    70%      3.501359e-04
    80%      4.001725e-04
    90%      4.500082e-04
    99%      4.950581e-04
    max      5.715329e-01
    Name: z_diff, dtype: float64
    Relative z diff stats:
    count    1.152128e+06
    mean     1.396460e-03
    std      8.912203e-03
    min      0.000000e+00
    10%      5.887793e-05
    20%      1.181598e-04
    30%      1.781700e-04
    40%      2.433816e-04
    50%      3.231510e-04
    60%      4.320990e-04
    70%      6.024769e-04
    80%      9.321848e-04
    90%      1.887206e-03
    99%      1.835998e-02
    max      6.995506e-01
    """
    df_peter = _prep_df_peter()
    df_peter= add_abs_zscore_to_df_peter(df_peter, suffix="_peter")
    df_martin = _prep_munged_df_martin()
    df_martin= _add_z_abs(df_martin)
    joined = df_martin.merge(df_peter, left_on="SNP", right_on="rsID")
    joined["z_diff"] = abs( abs(joined["Z_abs"]) - abs(joined["Z_abs_peter"])   )
    print("Absolute z diff stats:")
    print(joined["z_diff"].describe(percentiles=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]))
    print("Relative z diff stats:")
    joined = joined.loc[joined["Z_abs"]>0]
    joined["rel_z_diff"] = abs( (joined["Z_abs"] - joined["Z_abs_peter"] )/joined[["Z_abs","Z_abs_peter"]].max(axis=1))
    print(joined["rel_z_diff"].describe(percentiles=[0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]))
    print("yo")






if __name__ == "__main__":
    intersect_rsids()
    _get_value_counts()
    check_zscore_agreement()
    _get_allele_disagreement()

