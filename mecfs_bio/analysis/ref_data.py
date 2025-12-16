"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.db_snp.db_snp_build_37_as_parquet_unnest_ref import (
    PARQUET_DBSNP_37_UNNESTED,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW
            # DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED
            # LIU_ET_AL_2023_IBD_EUR_HARMONIZE
            # DB_SNP_VCF_FILE_BUILD_37        ,
            # PARQUET_DBSNP_37,
            PARQUET_DBSNP_37_UNNESTED
        ],
        incremental_save=True,
        must_rebuild_transitive=[PARQUET_DBSNP_37_UNNESTED],
    )


if __name__ == "__main__":
    run_initial_analysis()
