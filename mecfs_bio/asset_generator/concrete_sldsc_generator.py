"""
Asset generator for applying S-LDSC using standard reference datasets.
"""

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.custom.roadmap_cell_type_categorization import (
    ROADMAP_CELL_TYPE_CATEGORIES_FOR_LDSC,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_cahoy_ld_scores_extracted import (
    PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_corces_atac_ld_scores_extracted import (
    PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_EXTRACTED,
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
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.from_papers.finucane_2018_immgen_categories import (
    FICUANE_2018_IMMGEN_CATEGORIES,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.str_lowercase_pipe import StrLowercasePipe
from mecfs_bio.build_system.task.pipes.str_replace_pipe import StrReplacePipe
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe
from mecfs_bio.build_system.task_generator.sldsc_task_generator import (
    CellOrTissueLabelRecord,
    PartitionedLDScoresRecord,
    SLDSCTaskGenerator,
)


def standard_sldsc_task_generator(
    sumstats_task: Task,
    base_name: str,
    pre_pipe: DataProcessingPipe = IdentityPipe(),
) -> SLDSCTaskGenerator:
    """
    Standardized task generator for applying S-LSDC to GWAS summary statistics

    as described in
    "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types." Nature genetics 50.4 (2018): 621-629.

    using the standard partitioning datasets provided by the authors of that paper

    """
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
                cell_or_tissue_labels_task=CellOrTissueLabelRecord(
                    FICUANE_2018_IMMGEN_CATEGORIES,
                    pipe_left=CompositePipe(
                        [
                            StrLowercasePipe(
                                target_column="Name",
                                new_column_name="lower case name",
                            ),
                        ],
                    ),
                    pipe_right=CompositePipe(
                        [
                            StrLowercasePipe(
                                target_column="Cell Type",
                                new_column_name="lower case name",
                            )
                        ]
                    ),
                    left_join_on="lower case name",
                    right_join_on="lower case name",
                ),
                entry_name="immgen",
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="Corces_ATAC.ldcts",
                cell_or_tissue_labels_task=None,
                entry_name="corces_atac",
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="Multi_tissue_chromatin.ldcts",
                entry_name="multi_tissue_chromatin",
                cell_or_tissue_labels_task=CellOrTissueLabelRecord(
                    ROADMAP_CELL_TYPE_CATEGORIES_FOR_LDSC,
                    pipe_left=CompositePipe(
                        [
                            SplitExactColPipe(
                                col_to_split="Name",
                                split_by="__",
                                new_col_names=("Cell", "Epigenetic_Assay"),
                            ),
                            StrReplacePipe(
                                target_column="Cell",
                                new_column_name="Cell",
                                replace_what="_ENTEX",
                                replace_with="",
                            ),
                            StrLowercasePipe(
                                target_column="Cell",
                                new_column_name="Cell",
                            ),
                        ],
                    ),
                    pipe_right=CompositePipe(
                        [
                            StrLowercasePipe(
                                target_column="Cell type", new_column_name="Cell type"
                            )
                        ]
                    ),
                    left_join_on="Cell",
                    right_join_on="Cell type",
                ),
            ),
            PartitionedLDScoresRecord(
                ref_ld_chr_cts_task=PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED,
                ref_ld_chr_cts_filename="Multi_tissue_gene_expr.ldcts",
                cell_or_tissue_labels_task=CellOrTissueLabelRecord(
                    FICUANE_2018_FRANKE_GTEX_CATEGORIES,
                    pipe_left=CompositePipe(
                        [
                            StrLowercasePipe(
                                target_column="Name", new_column_name="Name"
                            ),
                            StrReplacePipe(
                                target_column="Name",
                                new_column_name="Name",
                                replace_what=" ",
                                replace_with="_",
                            ),
                        ]
                    ),
                    pipe_right=CompositePipe(
                        [
                            StrLowercasePipe(
                                target_column="Tissue_Or_Cell",
                                new_column_name="Tissue_Or_Cell",
                            ),
                            StrReplacePipe(
                                target_column="Tissue_Or_Cell",
                                new_column_name="Tissue_Or_Cell",
                                replace_what=" ",
                                replace_with="_",
                            ),
                        ]
                    ),
                ),
                entry_name="multi_tissue_gene_expression",
            ),
        ],
        multiple_testing_alpha=0.01,
        multiple_testing_method="fdr_bh",
        pre_pipe=pre_pipe,
    )
    return result
