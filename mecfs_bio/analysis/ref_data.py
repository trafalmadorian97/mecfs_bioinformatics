"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.processed.ukbb_rabgap1l_sumstats import (
    UKBBPPP_RABGAP1l_SUMSTATS_37_HARMONIZED,
)


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # UKBB_PPP_RABGAP1L,
            # UKBB_PPP_RABGAP1L_UNTAR,
            # STACK_UKBBPPP_RABGAP1L
            UKBBPPP_RABGAP1l_SUMSTATS_37_HARMONIZED
            # UNIPROT_LOOKUP
            # GENE_THESAURUS
            # DUNCAN_ET_AL_2025_ST1_EXTRACTED
        ],
        incremental_save=True,
        # must_rebuild_transitive=[SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED],
        # must_rebuild_transitive=[PARQUET_DBSNP_37_UNNESTED],
        # must_rebuild_transitive=[DB_SNP150_ANNOVAR_PROC]
    )


if __name__ == "__main__":
    run_initial_analysis()
