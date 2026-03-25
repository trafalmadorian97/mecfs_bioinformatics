from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    ManhattanPlotSettings,
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.raw.raw_xue_whole_brainstem import (
    XUE_WHOLE_BRAINSTEM_VOLUME_RAW,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)

XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        raw_gwas_data_task=XUE_WHOLE_BRAINSTEM_VOLUME_RAW,
        base_name="xue_whole_brainstem",
        fmt=GWASLabColumnSpecifiers(
            chrom="CHR",
            rsid="SNP",
            pos="POS",
            ea="A1",
            nea="A2",
            eaf="AF1",
            beta="BETA",
            se="SE",
            p="P",
        ),
        sample_size=83073,
        include_hba_magma_tasks=True,
        include_independent_cluster_plot_in_hba=True,
        manhattan_settings=ManhattanPlotSettings(),
        phenotype_info_for_ldsc=QuantPhenotype(),
        gget_settings=None,
    )
)
