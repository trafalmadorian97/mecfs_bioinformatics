"""
Genome-reference harmonization of the Liu et al. IBD sumstats (lifted over to build 37),
against the UCSC hg19 genome and the 1000 Genomes EUR panel.

- Orient every variant so that NEA is the plus-strand reference allele.
- Resolve palindromic SNVs and ambiguous indels by panel frequency when the table's
  orientation is not fully reference-consistent, and drop them when the evidence is not
  decisive.
"""

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats import (
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

LIU_ET_AL_2023_IBD_EUR_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
        asset_id="liu_et_al_2023_ibd_eur_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

LIU_ET_AL_2023_IBD_EUR_HARMONIZE = GenomeReferenceHarmonizationTask.create(
    asset_id="liu_et_al_2023_ibd_eur_genome_reference_harmonized",
    sumstats_task=LIU_ET_AL_2023_IBD_EUR_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
