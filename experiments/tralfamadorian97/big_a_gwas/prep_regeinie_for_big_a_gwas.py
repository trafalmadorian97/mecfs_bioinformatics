import polars as pl

from src.data_preprocessing_scripts.file_path_constants import (
    GWAS_1_AGWAS_PREP_PATH, GWAS_1_QC_APPLIED_FILE_PATH, GWAS_1_AGWAS_PREP_PATH_GZIPPED
)
from src.data_processing.preprocess_for_big_agwas import prep_for_big_agwas, BIG_A_GWAS_PVALUE_COL
from src.data_processing.decode_me_constants import DECODE_ME_SNP_COL
from src.data_processing_util.gzip_util import apply_gzip


def go():
    print("transforming")
    prep_for_big_agwas( GWAS_1_QC_APPLIED_FILE_PATH, GWAS_1_AGWAS_PREP_PATH)
    print("done transforming.")
    df = pl.scan_csv(GWAS_1_AGWAS_PREP_PATH, separator=" ")
    sig_genes = df.bottom_k(k=100, by=BIG_A_GWAS_PVALUE_COL).collect()
    print("significant genes in data before submission to big agwas")
    print(sig_genes)
    print("columns:")
    print(list(pl.scan_csv(GWAS_1_AGWAS_PREP_PATH, separator=" ").schema.keys()))
    print("gzipping")
    apply_gzip(GWAS_1_AGWAS_PREP_PATH, GWAS_1_AGWAS_PREP_PATH_GZIPPED)
    print("done gzipping")


if __name__ == "__main__":
    go()
