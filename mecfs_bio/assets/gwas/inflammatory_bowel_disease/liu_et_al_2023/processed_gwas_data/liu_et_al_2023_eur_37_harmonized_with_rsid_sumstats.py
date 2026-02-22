"""
Task to create a Gwaslab Sumstats objects from the Liu et al inflammatory bowel disease GWAS data with rsids assigned via the annovar version of dbSNP150
"""

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar import (
    LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.constants import gwaslab_constants

LIU_ET_AL_SUMSTATS_WITH_RSID_FROM_SNP150 = GWASLabCreateSumstatsTask(
    df_source_task=LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
    asset_id=AssetId("liu_et_al_2023_ibd_sumstats_with_rsids_via_snp150"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
    fmt=GWASLabColumnSpecifiers(
        nea=gwaslab_constants.GWASLAB_NON_EFFECT_ALLELE_COL,
        ea=gwaslab_constants.GWASLAB_EFFECT_ALLELE_COL,
        chrom=gwaslab_constants.GWASLAB_CHROM_COL,
        pos=gwaslab_constants.GWASLAB_POS_COL,
        eaf=gwaslab_constants.GWASLAB_EFFECT_ALLELE_FREQ_COL,
        beta=gwaslab_constants.GWASLAB_BETA_COL,
        se=gwaslab_constants.GWASLAB_SE_COL,
        p=gwaslab_constants.GWASLAB_P_COL,
        rsid="rsid",
        OR=None,
        info=None,
        snpid=gwaslab_constants.GWASLAB_SNPID_COL,
    ),
)
