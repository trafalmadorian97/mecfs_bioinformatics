import pandas as pd
from pathlib import Path


def process_names(path_to_mt_file:Path):
    df = pd.read_csv(path_to_mt_file)
    names_series = df["Name"]
    split_result = names_series.str.split("__", expand=True)
    split_result.columns= ["cell", "histone_mark"]
    split_result["cell"] = split_result["cell"].str.replace("_ENTEX","").str.title()
    cell_types = split_result["cell"].unique()
    ct = pd.Series(cell_types)
    ct.to_csv("cell_type_series.csv", index=False)
    import pdb; pdb.set_trace()
    print("yo")


if __name__ == "__main__":
    process_names(Path("assets/base_asset_store/gwas/ldl/million_veterans/analysis/million_veterans_eur_ldl_multi_tissue_chromatin_cell_analysis_multiple_testing_correction.csv"))