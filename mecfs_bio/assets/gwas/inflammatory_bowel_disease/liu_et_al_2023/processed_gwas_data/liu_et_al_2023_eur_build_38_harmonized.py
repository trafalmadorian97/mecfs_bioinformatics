"""
Genome-reference harmonization of the Liu et al. 2023 IBD European (NFE) summary statistics in
their native GRCh38 coordinates (no liftover), against the UCSC hg38 genome and the 1000 Genomes EUR (30x, GRCh38) allele-frequency panel.

The METAL-style alleles carry no reference orientation, so harmonization runs the untrusted
path: SNVs are oriented by the genome, palindromic SNVs and ambiguous indels by panel
frequency, or dropped.
"""

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_only import (
    LIU_ET_AL_2023_IBD_META_EUR_ONLY,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

LIU_ET_AL_2023_IBD_EUR_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=LIU_ET_AL_2023_IBD_META_EUR_ONLY,
    target_asset_id=AssetId("liu_et_al_2023_ibd_eur_build_38_sumstats"),
    basic_check=True,
    genome_build="38",
    fmt=GWASLabColumnSpecifiers(
        snpid="MarkerName",
        nea="Allele1",
        ea="Allele2",
        chrom="CHR",
        pos="BP",
        eaf="AF_NFE",
        beta="BETA_NFE",
        se="SE_NFE",
        p="P_NFE",
    ),
)

LIU_ET_AL_2023_IBD_EUR_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=LIU_ET_AL_2023_IBD_EUR_BUILD_38_SUMSTATS,
        asset_id="liu_et_al_2023_ibd_eur_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

LIU_ET_AL_2023_IBD_EUR_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="liu_et_al_2023_ibd_eur_build_38_genome_reference_harmonized",
    sumstats_task=LIU_ET_AL_2023_IBD_EUR_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
    panel_ancestry="eur",
)
