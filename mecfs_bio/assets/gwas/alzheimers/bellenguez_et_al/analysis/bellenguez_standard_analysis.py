from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import \
    concrete_standard_analysis_generator_assume_already_has_rsid
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.processed.extracted_bellinguez_data import \
    BELLINGUEZ_ET_AL_ALZHIEMERS_EXTRACTED
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.raw.raw_bellenguez_data import BELLENGUEZ_ET_AL_ALZHIEMERS_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabColumnSpecifiers
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe

BELLENGUEZ_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="bellenguez_et_al_alz",
    raw_gwas_data_task=BELLINGUEZ_ET_AL_ALZHIEMERS_EXTRACTED,
    sample_size=301478, # from gwas caltalog
    include_hba_magma_tasks=True,
    fmt=GWASLabColumnSpecifiers(
        rsid="hm_rsid",
        snpid="hm_variant_id",
        chrom="hm_chrom",
        pos="hm_pos",
        ea="hm_effect_allele",
        nea="hm_other_allele",
        OR="hm_odds_ratio",
        se="standard_error",
        p="p_value",
        info=None,
        eaf="hm_effect_allele_frequency",
        beta="hm_beta",
        ncase="n_cas",
        ncontrol="n_con",
    ),
    pre_pipe=DropNullsPipe(),
    hba_plot_settings=PlotSettings("plotly_white")
)