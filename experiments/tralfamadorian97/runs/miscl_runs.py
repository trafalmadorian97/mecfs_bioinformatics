from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_curated_gene_set_analysis import \
    DECODE_ME_CURATED_GENE_SET_ANALYSIS
from mecfs_bio.assets.gwas.multi_trait.genomic_sem.mecf_pain_common_factor import MECFS_PAIN_COMMON_FACTOR
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import \
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_gene_sets import \
    BENTHAM_2015_GENE_SET_ANALYSIS, BENTHAM_2015_GENE_SET_ANALYSIS_FROM_GENE_ANALYSIS
from mecfs_bio.assets.reference_data.gene_set_data.for_magma.from_gsea_msigdb.processed.full_msigdb_parquet_from_sqlite import MSIGDB_GENE_SETS_PARQUET_FROM_SQLLITE
from mecfs_bio.assets.reference_data.genomic_sem_reference.genomes_1k import GENOMES1K_REFERENCE_FOR_GENOMIC_SEM

from mecfs_bio.assets.reference_data.magma_specificity_matrices.processed.curated_potential_mecfs_gene_sets_specificity_matrix import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX
from mecfs_bio.assets.reference_data.magma_specificity_matrices.processed.curated_potential_mecfs_gene_sets_specificity_matrix_reduced import \
    CURATED_POTENTIAL_MECFS_GENE_SETS_SPECIFICITY_MATRIX_REDUCED
from mecfs_bio.assets.reference_data.genomic_sem_reference.hapmap3_snpist import HAPMAP3_SNPLIST_FOR_GENOMIC_SEM
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
        [
            MECFS_PAIN_COMMON_FACTOR
            # GENOMES1K_REFERENCE_FOR_GENOMIC_SEM
            # HAPMAP3_SNPLIST_FOR_GENOMIC_SEM
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
        ],
        settings=TopologicalSchedulerSettings(
            print_progress=False
        )
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_miscl_analysis()
