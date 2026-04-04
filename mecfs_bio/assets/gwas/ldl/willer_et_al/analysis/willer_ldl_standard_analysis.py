from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import \
    concrete_standard_analysis_generator_assume_already_has_rsid
from mecfs_bio.assets.gwas.ldl.willer_et_al.raw.raw_willer_ldl_data import WILLER_LDL_EUR_DATA_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import QuantPhenotype

WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS=concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="willer_et_al_2023_ldl_eur",
    raw_gwas_data_task=WILLER_LDL_EUR_DATA_RAW,
    fmt=GWASLabColumnSpecifiers(
        rsid="hm_rsid",
        chrom="hm_chrom",
        pos="hm_pos",
        nea="hm_other_allele",
        ea="hm_effect_allele",
        beta="hm_beta",
        se="standard_error",
        p="p_value",
        eaf="eaf_ref"
    ),
    sample_size=89138,
   phenotype_info_for_ldsc=QuantPhenotype(),

)