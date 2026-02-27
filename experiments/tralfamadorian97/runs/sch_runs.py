from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.sch_clustermap_plot import \
    SCH_PGC_2022_MAGMA_CLUSTER_MAP


def run_initial_schizophrenia_analysis():
    DEFAULT_RUNNER.run(
        [
            SCH_PGC_2022_MAGMA_CLUSTER_MAP
            # SCH_PGC_2022_MAGMA_CLUSTER_MAP
            # SCH_PGC_2022_MAGMA_CLUSTER_MAP_PLOTLY
            # SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.gene_analysis_task,
            # SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.filtered_gene_analysis_task,
            # GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA
        ],
        must_rebuild_transitive=[
            SCH_PGC_2022_MAGMA_CLUSTER_MAP
            # SCH_PGC_2022_MAGMA_CLUSTER_MAP
            # SCH_PGC_2022_MAGMA_CLUSTER_MAP_PLOTLY
        ],
    )


if __name__ == "__main__":
    run_initial_schizophrenia_analysis()
