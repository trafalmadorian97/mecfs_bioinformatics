from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.processed.full_msigdb_parquet_from_sqlite import MSIGDB_GENE_SETS_PARQUET_FROM_SQLLITE

from mecfs_bio.assets.reference_data.magma_specificity_matrices.processed.curated_potential_mecfs_gene_sets_specificity_matrix import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX
from mecfs_bio.assets.reference_data.magma_specificity_matrices.processed.curated_potential_mecfs_gene_sets_specificity_matrix_reduced import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX_REDUCED
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
        [
            CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX_REDUCED
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
        ],
        settings=TopologicalSchedulerSettings(
            print_progress=False
        )
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
