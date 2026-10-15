"""
Genome-reference harmonization of the UK Biobank ankylosing spondylitis sumstats (lifted
over to build 37), against the UCSC hg19 genome and the gnomAD v2.1.1 genomes panel, so
that NEA is the plus-strand reference allele. Ambiguous variants are resolved with the
gnomAD north-western European (nfe_nwe) allele frequencies, the subgroup covering UK
cohorts.

Duplicated variant keys, likely a consequence of liftover, are reduced to the row with
the smallest SE before harmonization.
"""

from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_sumstats import (
    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

UKBB_ANK_SPOND_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS,
        asset_id="uk_biobank_ank_spond_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

UKBB_ANK_SPOND_GENOME_REFERENCE_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="uk_biobank_ank_spond_genome_reference_harmonized",
    sumstats_task=UKBB_ANK_SPOND_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    panel_task=GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
    panel_ancestry="nfe_nwe",
    pipe=UniquePipe(
        by=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],
        keep="first",
        order_by=[GWASLAB_SE_COL],
    ),
)
