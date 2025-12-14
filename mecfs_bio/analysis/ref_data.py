"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.db_snp.db_snp_vcf_file_with_index import (
    DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_analysis()
