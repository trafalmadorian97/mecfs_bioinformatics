from mecfs_bio.assets.gwas.ldl.million_veterans.processed_gwas_data.million_veterans_ldl_eur_magma_task_generator import \
    MILLION_VETERANS_EUR_LDL_MAGMA_TASKS
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_gene_expr_ld_score_extracted import \
    PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_regression_weights_extracted import \
    PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_extracted import \
    BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_ldsc import CellAnalysisByLDSCTask

MILLION_VETERAN_LDL_BASE_CELL_ANALYSIS_BY_LDSC = CellAnalysisByLDSCTask.create(
   "million_veterans_ld_base_cell_analysis_by_ldsc" ,
    source_sumstats_task=MILLION_VETERANS_EUR_LDL_MAGMA_TASKS.sumstats_task,
    ref_ld_chr_cts_task=PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED,
    ref_ld_chr_cts_filename="Multi_tissue_gene_expr.ldcts",
    ref_ld_chr_task=BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
    ref_ld_chr_inner_dirname="baseline_v1.2/baseline.@",
    w_ld_chr_task=PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
    w_ld_chr_inner_dirname="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@"
    
)