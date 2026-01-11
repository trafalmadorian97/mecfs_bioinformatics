import pandas as pd
from pathlib import Path


def load_immgen_cell_types(path: Path):
    table = pd.read_csv(path)
    split_result = table["Name"].str.split(".",expand=True)
    cell_types = split_result[0].unique()
    print(cell_types)
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == "__main__":
    load_immgen_cell_types(Path("assets/base_asset_store/gwas/ldl/million_veterans/analysis/million_veterans_eur_ldl_immgen_cell_analysis_multiple_testing_correction.csv"))