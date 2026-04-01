from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_only import (
    LIU_ET_AL_2023_IBD_META_EUR_ONLY,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)

LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=LIU_ET_AL_2023_IBD_META_EUR_ONLY,
    target_asset_id=AssetId("liu_et_al_2023_ibd_sumstats_liftover_to_37"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
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
        rsid=None,
        OR=None,
        info=None,
    ),
)
