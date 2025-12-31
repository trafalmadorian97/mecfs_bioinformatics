from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import \
    concrete_standard_analysis_generator_assume_already_has_rsid
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.raw.raw_dbp_data import KEATON_ET_AL_DBP_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers

KEATON_DBP_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="keaton_dbp",
    raw_gwas_data_task=KEATON_ET_AL_DBP_RAW,
    fmt = GWASLabColumnSpecifiers(
    rsid="rsid",
    chrom="chromosome",
    pos="base_pair_location",
    ea="effect_allele",
    nea="other_allele",
    beta="beta",
    eaf="effect_allele_frequency",
    p="p_value",
    # OR="odds_ratio",
    se="standard_error",
    n="n",
    # ncase="num_cases",
    # ncontrol="num_controls",
    # info="r2",
    snpid="variant_id",
    # snpid=None,
    # or_95l="ci_lower",
    # or_95u="ci_upper",
        OR=None,
        info=None
),
    sample_size=729882, # from summary statistics file
    include_master_gene_lists=False
)