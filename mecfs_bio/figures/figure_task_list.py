"""
This file contains a list of figures to be exported for use in documentation
"""

from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.analysis.xue_whole_brainstem_standard_analysis import (
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.c_reactive_protein.said_et_al.analysis.said_crp_standard_analysis import (
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.ebv_dna.nyeo_et_al.analysis.ebv_dna_standard_analysis import (
    EBV_DNA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.height.yengo_2022.analysis.yengo_standard_analysis import (
    YENGO_HEIGHT_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.human_herpesvirus_7_dna.kamitaki_et_al_2025.analysis.kamitaki_et_al_2025_standard_analysis import (
    KAMITAKI_ET_AL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.mv_ldl_heritability_task import (
    MV_LDL_LDSC_RESULTS_MARKDOWN,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_standard_analysis import (
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc import (
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_manhattan import (
    DECODE_ME_GWAS_1_MANHATTAN_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_curated_gene_set_analysis import (
    DECODE_ME_CURATED_GENE_SET_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_magma_gene_plot import (
    DECODE_ME_MAGMA_GENE_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.mixer.decode_me_univariate_mixer import (
    DECODE_ME_UNIVARIATE_MIXER,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.analysis.million_veterans_cfs_standard_analysis import (
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies import (
    CFS_CT_LDSC_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.analysis.neale_lab_cfs_standard_analysis import (
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.migraine.million_veterans.analysis.million_veterans_migraine_standard_analysis import (
    MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_plot import (
    CT_LDSC_INITIAL_PLOT,
)
from mecfs_bio.assets.gwas.multi_trait.lcv.mi_lcv_analysis import MI_LCV_TASK_GROUP
from mecfs_bio.assets.gwas.multi_trait.polygenic_overlap.bivariate_mixer.mecfs_pain_bivariate_mixer import (
    MECFS_PAIN_BIVARIATE_MIXER,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.mixer.johnston_et_al_univariate_mixer import (
    JOHNSTON_ET_AL_UNIVARIATE_MIXER,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.triglycerides.willer_et_al.analysis.triglycide_standard_analysis import (
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.base_task import Task

MULTI_TISSUE_CHROMATIN_REF = "multi_tissue_chromatin"
MULTI_TISSUE_GENE_EXPRESSION_REF = "multi_tissue_gene_expression"

ALL_FIGURE_TASKS: list[Task] = [
    DECODE_ME_GWAS_1_MANHATTAN_PLOT,
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD,
    DECODE_ME_MAGMA_GENE_PLOT,
    CT_LDSC_INITIAL_PLOT,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.extracted_plot_task,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task_unwrap,
    BENTHAM_LUPUS_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_GENE_EXPRESSION_REF
    ].plot_task_unwrap,
    EBV_DNA_STANDARD_ANALYSIS.tasks.hba_magma_tasks_unwrap.extracted_plot_task,
    EBV_DNA_STANDARD_ANALYSIS.tasks.magma_tasks.inner.bar_plot_task,
    EBV_DNA_STANDARD_ANALYSIS.tasks.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task_unwrap,
    DECODE_ME_UNIVARIATE_MIXER.power_plot_task,
    DECODE_ME_UNIVARIATE_MIXER.qq_plot_task,
    DECODE_ME_UNIVARIATE_MIXER.result_markdown_table_task,
    DECODE_ME_CURATED_GENE_SET_ANALYSIS.bar_plot_task_full,
    JOHNSTON_ET_AL_UNIVARIATE_MIXER.power_plot_task,
    JOHNSTON_ET_AL_UNIVARIATE_MIXER.qq_plot_task,
    JOHNSTON_ET_AL_UNIVARIATE_MIXER.result_markdown_table_task,
    KAMITAKI_ET_AL_STANDARD_ANALYSIS.tasks.hba_magma_tasks_unwrap.extracted_plot_task,
    KAMITAKI_ET_AL_STANDARD_ANALYSIS.tasks.hba_magma_tasks_unwrap.independent_clusters_markdown_task_unwrap,
    KAMITAKI_ET_AL_STANDARD_ANALYSIS.tasks.magma_tasks.inner.bar_plot_task,
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.extracted_plot_task,
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.independent_clusters_markdown_task_unwrap,
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.magma_independent_cluster_plot_unwrap,
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task_unwrap,
    MECFS_PAIN_BIVARIATE_MIXER.result_table_markdown_task,
    MV_LDL_LDSC_RESULTS_MARKDOWN,
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_GENE_EXPRESSION_REF
    ].plot_task_unwrap,
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_CHROMATIN_REF
    ].plot_task_unwrap,
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.sldsc_tasks.partitioned_tasks[
        MULTI_TISSUE_GENE_EXPRESSION_REF
    ].plot_task_unwrap,
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task,
    MI_LCV_TASK_GROUP.downstream_trait_tables["MI"],
    YENGO_HEIGHT_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
    MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS.heritability_markdown_task_unwrap,
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.heritability_markdown_task_unwrap,
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.magma_tasks.inner.bar_plot_task,
    CFS_CT_LDSC_ASSET_GENERATOR.aggregation_markdown_task,
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP.heritability_markdown_task_unwrap,
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP.magma_tasks.inner.bar_plot_task,
]
