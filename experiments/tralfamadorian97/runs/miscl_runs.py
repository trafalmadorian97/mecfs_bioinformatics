from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.height.yengo_2022.analysis.yengo_standard_analysis import YENGO_HEIGHT_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.ldl.multistudy.genetic_correlation.ct_ldsc.ct_ldsc_ldl import CT_LDSC_LDL
from mecfs_bio.assets.gwas.me_cfs.astra_zenica_phewas_gene_level.raw.get_mecfs_az_phewas import MECFS_AZ_PHEWAS
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc import DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC, \
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_magma_gene_plot import DECODE_ME_MAGMA_GENE_PLOT
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_magma_gene_plot_with_window import \
    DECODE_ME_MAGMA_GENE_PLOT_WITH_WINDOW
from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies import \
    CFS_CT_LDSC_ASSET_GENERATOR
from mecfs_bio.assets.gwas.me_cfs.multistudy.analysis.genetic_correlation.ct_ldsc.ct_ldsc_mecfs_studies_plot import \
    CT_LDSC_CFS_CORR_PLOT
from mecfs_bio.build_system.scheduler.topological_scheduler import TopologicalSchedulerSettings


def run_miscl_analysis():
    DEFAULT_RUNNER.run(

            # WILLER_LDL_EUR_DATA_RAW
            # MV_LDL_LDSC_RESULTS_MARKDOWN
        # MV_LDL_HERITABILITY_TASK
        #     MI_LDL_WILLER_CORRELATION.terminal_tasks()
    # [
        # LDL_MI_LCV_ANALYSIS,
        # [MI_LCV_PLOT]
        # [

            # YENGO_HEIGHT_STANDARD_ANALYSIS.heritability_markdown_task_unwrap
            # MILLION_VETERAN_MIGRAINE_EUR_DATA_RAW,
            # MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS.get_terminal_tasks()
        # BENTHAM_2015_GENE_SET_ANALYSIS.terminal_tasks()+
        # BENTHAM_2015_GENE_SET_ANALYSIS_FROM_GENE_ANALYSIS.terminal_tasks()+
        # JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.gene_set_analysis_tasks.terminal_tasks()+

        # DECODE_ME_CURATED_GENE_SET_ANALYSIS.terminal_tasks()+
        # MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP.get_terminal_tasks()+
        # CFS_CT_LDSC_ASSET_GENERATOR.terminal_tasks()+
        # CT_LDSC_LDL.terminal_tasks()+
        [
            DECODE_ME_MAGMA_GENE_PLOT,
            DECODE_ME_MAGMA_GENE_PLOT_WITH_WINDOW,
            # YENGO_HEIGHT_STANDARD_ANALYSIS.magma_tasks.inner.bar_plot_task
            # MECFS_AZ_PHEWAS
            # CT_LDSC_CFS_CORR_PLOT
            # MILLION_VETERANS_CFS_RAW,
            # DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC
            # DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD
            # CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX_REDUCED
            # MSIGDB_GENE_SETS_PARQUET_FROM_SQLLITE,
            # MSIGDB_SQLLITE_EXTRACTED
            # MSIGDB_GENE_SETS_PARQUET_UNPACKED
            # MSIGDB_JSON_GENE_SETS_PARQUET
            # MSIGDB_JSON_GENE_SETS
        ]

        # ]
         # YENGO_STANDARD_ANALYSIS.get_terminal_tasks()
        # WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS.get_terminal_tasks()
        # WILLER_TG_EUR_DATA_RAW
     # MECFS_PAIN_LCV_ANALYSIS
     # ]
            ,
        # +WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_DATA_RAW]+MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task],
        # [MI_EUR_MANHATTAN],

        incremental_save=True,
        must_rebuild_transitive=[
            # CFS_CT_LDSC_ASSET_GENERATOR.aggregation_markdown_task
        ],
        settings=TopologicalSchedulerSettings(
            print_progress=False
        )
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
