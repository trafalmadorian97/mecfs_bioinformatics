"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW
            # DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR
            # MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED
            # LIU_ET_AL_2023_IBD_EUR_HARMONIZE
            # DB_SNP_VCF_FILE_BUILD_37        ,
            # PARQUET_DBSNP_37,
            # PARQUET_DBSNP_37_UNNESTED,
            # PARQUET_DBSNP_37_UNNESTED_SORTED
            # DB_SNP150_ANNOVAR_PROC,
            # PARQUET_DBSNP150_37_ANNOVAR_PROC,
            # PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME
            # PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE
        ],
        incremental_save=True,
        # must_rebuild_transitive=[PARQUET_DBSNP_37_UNNESTED],
        # must_rebuild_transitive=[DB_SNP150_ANNOVAR_PROC]
    )


if __name__ == "__main__":
    run_initial_analysis()
