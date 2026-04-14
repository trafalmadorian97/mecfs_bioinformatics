from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.height.yengo_2022.raw.yengo_height_raw import YENGO_ET_AL_RAW
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    QuantPhenotype,
)
from mecfs_bio.build_system.task.magma.magma_plot_brain_atlas_result_with_stepwise_labels import (
    HBAIndepPlotOptions,
)
from mecfs_bio.build_system.task.magma.plot_magma_brain_atlas_result import PlotSettings
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import (
    FilterRowsByMinInCol,
)

YENGO_HEIGHT_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="yengo_et_al_height",
        raw_gwas_data_task=YENGO_ET_AL_RAW,
        sample_size=int(1.587709e06),  # sample size varies by SNP.  This is the median
        fmt=GWASLabColumnSpecifiers(
            rsid="RSID",
            chrom="CHR",
            pos="POS",
            ea="EFFECT_ALLELE",
            nea="OTHER_ALLELE",
            eaf="EFFECT_ALLELE_FREQ",
            beta="BETA",
            p="P",
            se="SE",
        ),
        include_hba_magma_tasks=True,
        hba_plot_settings=PlotSettings("plotly_white"),
        include_independent_cluster_plot_in_hba=True,
        hba_indep_plot_options=HBAIndepPlotOptions(annotation_text_size=11),
        pre_pipe=FilterRowsByMinInCol(
            col="N",
            min_value=1.032694e06,
        ),
        phenotype_info_for_ldsc=QuantPhenotype(),
    )
)
