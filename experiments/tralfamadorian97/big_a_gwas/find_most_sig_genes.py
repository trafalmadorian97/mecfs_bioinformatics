import polars as pl

from src.data_preprocessing_scripts.file_path_constants import GWAS_1_AGWAS_PREP_PATH
from src.data_processing.preprocess_for_big_agwas import BIG_A_GWAS_PVALUE_COL


def most_sig_genes():
    df = pl.scan_csv(GWAS_1_AGWAS_PREP_PATH, separator=" ")
    sig_genes = df.bottom_k(k=100, by=BIG_A_GWAS_PVALUE_COL).collect()
    print(sig_genes)


if __name__ == "__main__":
    most_sig_genes()
