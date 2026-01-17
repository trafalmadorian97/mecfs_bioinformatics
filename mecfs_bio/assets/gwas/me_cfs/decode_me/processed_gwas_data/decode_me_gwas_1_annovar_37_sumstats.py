from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.constants import gwaslab_constants

DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    asset_id=AssetId("decode_me_gwas_1_37_sumstats_rsids_from_annovar"),
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
        p=None,
        chi_sq=gwaslab_constants.GWASLAB_CHISQ_COL,
        mlog10p=gwaslab_constants.GWASLAB_MLOG10P_COL,
        rsid="rsid",
        OR=None,
        info=None,
        snpid=gwaslab_constants.GWASLAB_SNPID_COL,
        n=gwaslab_constants.GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
)
