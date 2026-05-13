from mecfs_bio.asset_generator.magma_curated_gene_set_analysis_generator import (
    curated_gene_set_analysis_magma_tasks_from_gene_analysis,
    generate_curated_gene_set_analysis_magma_tasks,
)
from mecfs_bio.assets.gwas.systemic_lupus_erythematosus.bentham_et_al_2015.analysis_results.bentham_2015_standard_analysis import (
    BENTHAM_LUPUS_STANDARD_ANALYSIS,
)

BENTHAM_2015_GENE_SET_ANALYSIS = generate_curated_gene_set_analysis_magma_tasks(
    base_name="betham_2015",
    gwas_parquet_with_rsids_task=BENTHAM_LUPUS_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
    sample_size=14267,
)
BENTHAM_2015_GENE_SET_ANALYSIS_FROM_GENE_ANALYSIS = curated_gene_set_analysis_magma_tasks_from_gene_analysis(
    base_name="bentham_2015",
    gene_analysis_task=BENTHAM_LUPUS_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.magma_gene_analysis_task,
)
