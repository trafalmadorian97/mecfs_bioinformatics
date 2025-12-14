"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW
            # DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED
            LIU_ET_AL_2023_IBD_EUR_HARMONIZE
        ],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_analysis()
