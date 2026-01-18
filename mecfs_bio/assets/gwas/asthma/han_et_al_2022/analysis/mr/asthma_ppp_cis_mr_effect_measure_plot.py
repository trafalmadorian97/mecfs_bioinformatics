from mecfs_bio.assets.gwas.asthma.han_et_al_2022.analysis.mr.asthma_ppp_cis_mr_multiple_testing import (
    HAN_2022_ASTHMA_CIS_PPP_TSMR_MULTIPLE_TESTING,
)
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe
from mecfs_bio.build_system.task.plot_mr_effect_measure_task import (
    EffectMeasurePlotConfig,
    PlotMREffectMeasure,
)
from mecfs_bio.build_system.task.two_sample_mr_task import (
    TSM_OUTPUT_B_COL,
    TSM_OUTPUT_EXPOSURE_COL,
    TSM_OUTPUT_SE_COL,
)
from mecfs_bio.constants.sun_et_al_pqtl_constants import (
    SUN_ASSAY_TARGET,
    SUN_TARGET_UNIPROT,
)

HAN_CIS_PPP_MR_PLOT = PlotMREffectMeasure.create(
    asset_id="cis_ppp_mr_effect_measure_plot",
    source_df_task=HAN_2022_ASTHMA_CIS_PPP_TSMR_MULTIPLE_TESTING,
    config=EffectMeasurePlotConfig(
        y_label_col=SUN_ASSAY_TARGET,
        y_label="Protein",
        effect_size_col=TSM_OUTPUT_B_COL,
        effect_size_label="b",
        se_col=TSM_OUTPUT_SE_COL,
        ref_line_center=0,
        figsize=(14 * 1.3, 6 * 1.3),
        t_adjust=0.02,
    ),
    pre_pipe=SplitExactColPipe(
        TSM_OUTPUT_EXPOSURE_COL,
        split_by="_",
        new_col_names=(SUN_ASSAY_TARGET, SUN_TARGET_UNIPROT),
    ),
)
