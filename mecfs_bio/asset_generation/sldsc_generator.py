from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_cahoy_ld_scores_extracted import (
    PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_gtex_brain_ld_scores_extracted import (
    PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_immgen_ld_scores_extracted import (
    PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_chromatin_ld_scores_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_gene_expr_ld_score_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_regression_weights_extracted import (
    PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_extracted import (
    BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.from_papers.finucane_2018_franke_gtex_categories import (
    FICUANE_2018_FRANKE_GTEX_CATEGORIES,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task_generator.sldsc_task_generator import (
    PartitionedLDScoresRecord,
    SLDSCTaskGenerator,
)


def standard_sldsc_task_generator(
    sumstats_task: Task, base_name: str
) -> SLDSCTaskGenerator:
    result = SLDSCTaskGenerator.create(
        base_name=base_name,
        source_sumstats_task=sumstats_task,
        ref_ld_chr_task=BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
        ref_ld_chr_inner_dirname="baseline_v1.2/baseline.@",
        w_ld_chr_task=PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
        w_ld_chr_inner_dirname="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@",
        partitioned_entries=[
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="Cahoy.ldcts",
                cell_or_tissue_labels_task=None,
                entry_name="cahoy_cns",
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="GTEx_brain.ldcts",
                cell_or_tissue_labels_task=None,
                entry_name="gtex_brain",
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="ImmGen.ldcts",
                cell_or_tissue_labels_task=None,
                entry_name="immgen",
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="Multi_tissue_chromatin.ldcts",
                cell_or_tissue_labels_task=None,
                entry_name="multi_tissue_chromatin",
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="Multi_tissue_gene_expr.ldcts",
                cell_or_tissue_labels_task=FICUANE_2018_FRANKE_GTEX_CATEGORIES,
                entry_name="multi_tissue_gene_expression",
            ),
        ],
        multiple_testing_alpha=0.05,
        multiple_testing_method="fdr_bh",
    )
    return result
