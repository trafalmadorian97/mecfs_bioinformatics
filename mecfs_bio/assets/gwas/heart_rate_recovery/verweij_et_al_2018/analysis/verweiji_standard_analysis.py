"""
Task to perform standard analysis on heart rate recovery GWAS of Verweiji et al.
"""

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.heart_rate_recovery.verweij_et_al_2018.raw.verweij_raw_summary import (
    VERWEIJI_ET_AL_RAW_HARMONIZED_BUILD_37,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings

VERWEIJI_ET_AL_HRR_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="verweiji_et_al_hrr",
        raw_gwas_data_task=VERWEIJI_ET_AL_RAW_HARMONIZED_BUILD_37,
        sample_size=58_818,  # from abstract
        fmt=GWASLabColumnSpecifiers(
            rsid="variant_id",
            chrom="chromosome",
            pos="base_pair_location",
            ea="effect_allele",
            nea="other_allele",
            beta="beta",
            p="p_value",
            snpid="uniqid",
            OR=None,
            se="standard_error",
            info="info",
        ),
        include_hba_magma_tasks=True,
        hba_plot_settings=PlotSettings("plotly_white"),
        include_independent_cluster_plot_in_hba=True,
        hba_indep_plot_options=HBAIndepPlotOptions(annotation_text_size=11),
    )
)
