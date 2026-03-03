from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.raw.raw_mv_eur_ank_spond_data import \
    MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabCreateSumstatsTask, \
    GWASLabColumnSpecifiers

MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_SUMSTATS=GWASLabCreateSumstatsTask(
    df_source_task=MILLION_VETERAN_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW,
    asset_id="million_veterans_ank_spond_sumstats",
    fmt=GWASLabColumnSpecifiers(
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        OR="odds_ratio",
        eaf="effect_allele_frequency",
        p="p_value",
        rsid="rsid",
        or_95u="ci_upper",
        or_95l="ci_lower",
        ncase="num_cases",
        ncontrol="num_controls"
    ),
    liftover_to="19",
    genome_build="infer",
    basic_check=True
)