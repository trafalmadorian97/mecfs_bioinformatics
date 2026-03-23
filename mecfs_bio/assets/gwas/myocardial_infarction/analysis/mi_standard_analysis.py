from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    ManhattanPlotSettings,
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.myocardial_infarction.raw.raw_mi_data import (
    MILLION_VETERAN_MI_EUR_DATA_RAW,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe

MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS = concrete_standard_analysis_generator_assume_already_has_rsid(
    base_name="million_veterans_mi",
    raw_gwas_data_task=MILLION_VETERAN_MI_EUR_DATA_RAW,
    sample_size=432053,  # from GWAS catalog yaml file
    fmt=GWASLabColumnSpecifiers(
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        OR="odds_ratio",
        eaf="effect_allele_frequency",
        p="p_value",
        rsid="rsid",
    ),
    include_hba_magma_tasks=True,
    hba_plot_settings=PlotSettings("plotly_white"),
    include_independent_cluster_plot_in_hba=True,
    pre_sldsc_pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    include_master_gene_lists=False,
    gget_settings=None,
    manhattan_settings=ManhattanPlotSettings(),
    phenotype_info_for_ldsc=BinaryPhenotypeSampleInfo(
        sample_prevalence=0.09,  # 39,074/(39,074+392,979)
        estimated_population_prevalence=0.04,  # https://www.ncbi.nlm.nih.gov/books/NBK83160/
    ),
)
