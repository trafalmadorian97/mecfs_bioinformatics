from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.raw_gwas_data.bentham_2015_raw_build_37 import (
    BENTHAM_2015_RAW_BUILD_37,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings

BENTHAM_LUPUS_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        "bentham_2015_lupus",
        raw_gwas_data_task=BENTHAM_2015_RAW_BUILD_37,
        fmt=GWASLabColumnSpecifiers(
            rsid="rsid",
            snpid=None,
            chrom="chrom",
            pos="pos",
            ea="effect_allele",
            nea="other_allele",
            OR="OR",
            se="se",
            p="p",
            info=None,
            beta="beta",
        ),
        sample_size=14267,  # from gwas catalog metadata yaml file
        include_hba_magma_tasks=True,
        include_independent_cluster_plot_in_hba=True,
        hba_plot_settings=PlotSettings("plotly_white"),
    )
)
