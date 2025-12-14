"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.red_blood_cell_volume.million_veterans.raw_gwas_data.raw_red_blood_volume import (
    MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW
            # DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR
            MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED
        ],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_analysis()
