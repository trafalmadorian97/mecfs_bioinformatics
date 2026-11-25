"""
Prepared (parquet) forms of the standard S-LDSC reference data, for the vectorized
cell-type-specific S-LDSC Task: the six LDSC-SEG annotation sets of Finucane et al. 2018 and the
baseline v1.2 LD scores with HapMap3 (no MHC) regression weights.
"""

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_cahoy_ld_scores_extracted import (
    PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_corces_atac_ld_scores_extracted import (
    PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_gtex_brain_ld_scores_extracted import (
    PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_immgen_ld_scores_extracted import (
    PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_chromatin_ld_scores_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_multi_tissue_gene_expr_ld_score_extracted import (
    PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED_RD,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.partitioned_model_regression_weights_extracted import (
    PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.thousand_genome_baseline_ld_extracted import (
    BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
)
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PrepareLdctsAnnotationSetTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_sldsc_baseline_task import (
    PrepareSLDSCBaselineTask,
)

PREPARED_BASELINE_V1_2_HM3_WEIGHTS = PrepareSLDSCBaselineTask.create(
    asset_id="prepared_sldsc_baseline_v1_2_hm3_no_mhc_weights",
    baseline_directory_task=BASE_MODEL_PARTITIONED_LD_SCORES_EXTRACTED,
    baseline_pattern="baseline_v1.2/baseline.@",
    weights_directory_task=PARTITIONED_MODEL_REGRESSION_WEIGHTS_EXTRACTED,
    weights_pattern="1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.@",
)

PREPARED_CAHOY_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_cahoy_annotation_set",
    source_directory_task=PARTITIONED_MODEL_CAHOY_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Cahoy.ldcts",
)
PREPARED_GTEX_BRAIN_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_gtex_brain_annotation_set",
    source_directory_task=PARTITIONED_MODEL_GTEX_BRAIN_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="GTEx_brain.ldcts",
)
PREPARED_IMMGEN_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_immgen_annotation_set",
    source_directory_task=PARTITIONED_MODEL_IMMGEN_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="ImmGen.ldcts",
)
PREPARED_CORCES_ATAC_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_corces_atac_annotation_set",
    source_directory_task=PARTITIONED_MODEL_CORCES_ATAC_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Corces_ATAC.ldcts",
)
PREPARED_MULTI_TISSUE_CHROMATIN_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_multi_tissue_chromatin_annotation_set",
    source_directory_task=PARTITIONED_MODEL_MULTI_TISSUE_CHROMATIN_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Multi_tissue_chromatin.ldcts",
)
PREPARED_MULTI_TISSUE_GENE_EXPR_ANNOTATION_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prepared_sldsc_multi_tissue_gene_expr_annotation_set",
    source_directory_task=PARTITIONED_MODEL_MULTI_TISSUE_GENE_EXPR_LD_SCORES_EXTRACTED_RD,
    ldcts_filename="Multi_tissue_gene_expr.ldcts",
)
