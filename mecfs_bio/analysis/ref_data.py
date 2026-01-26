"""
Script to download some reference data.  Mainly for testing.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import GENE_THESAURUS
from mecfs_bio.assets.reference_data.uniprot.uniprot_lookup_table import UNIPROT_LOOKUP


def run_initial_analysis():
    DEFAULT_RUNNER.run(
        [
            # UNIPROT_LOOKUP
            GENE_THESAURUS
            # DUNCAN_ET_AL_2025_ST1_EXTRACTED
        ],
        incremental_save=True,
        # must_rebuild_transitive=[SUN_ET_AL_2023_COMBINED_PQTLS_EXTRACTED],
        # must_rebuild_transitive=[PARQUET_DBSNP_37_UNNESTED],
        # must_rebuild_transitive=[DB_SNP150_ANNOVAR_PROC]
    )


if __name__ == "__main__":
    run_initial_analysis()
