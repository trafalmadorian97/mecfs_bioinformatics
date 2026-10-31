"""
Genome-reference harmonization of the DecodeME GWAS-1 summary statistics (native GRCh38),
against the UCSC hg38 genome and the 1000 Genomes EUR (30x, GRCh38) allele-frequency panel.

- Orient every variant so that NEA is the plus-strand GRCh38 reference allele.
- Resolve palindromic SNVs and ambiguous indels by panel frequency when the table's
  orientation is not fully reference-consistent, and drop them when the evidence is not
  decisive.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_minimal_processing import (
    DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

DECODE_ME_GWAS_1_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
        asset_id="decode_me_gwas_1_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

DECODE_ME_GWAS_1_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="decode_me_gwas_1_build_38_genome_reference_harmonized",
    sumstats_task=DECODE_ME_GWAS_1_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
    panel_ancestry="eur",
)
