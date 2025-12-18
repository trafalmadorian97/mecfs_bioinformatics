"""
Script to analyze inflammatory bowl disease GWAS from Liu et al.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar import (
    LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
)


def run_ibd_analysis():
    """
    This script runs a basic analysis of the IBD GWAS from Liu et al.
    It includes:


    - Using MAGMA to combine the IBD GWAS summary statistics with tissue-specific expression data from GTEx to identify possible tissues involved in the IBD disease process.
    """
    DEFAULT_RUNNER.run(
        [
            # LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
            # LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
            #    LIU_ET_AL_2023_ASSIGN_RSID_VIA_DBSNP,
            # LIU_ET_AL_2023_ASSIGN_RSID_VIA_DBSNP_NON_HARMONIZED,
            LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR
        ],
        incremental_save=True,
        must_rebuild_transitive=[LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR],
    )


if __name__ == "__main__":
    run_ibd_analysis()
