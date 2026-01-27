"""
Generate a clustermap showing gene/tissue specificity for the most significant Schizophrenia genes as detected by MAGMA
"""

from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.standard_analysis_sc_pgc_2022 import (
    SCH_PGC_2022_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.assets.reference_data.rna_seq_data.gtex_v10_median_tissue_expression_rna_seq_prep_for_magma import (
    GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
)
from mecfs_bio.build_system.data_manipulation.xr_data.xr_gene_dataset_load import (
    GeneInfoSource,
    SpecificityMatrixSource,
)
from mecfs_bio.build_system.task.gene_tissue_expression_clustermap_task import (
    ExpressionMatrixClusterMapTaskV2,
)
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.magma_constants import MAGMA_GENE_COL

SCH_PGC_2022_MAGMA_CLUSTER_MAP = ExpressionMatrixClusterMapTaskV2.create_standard_gene_magma_heatmap(
    asset_id="sch_pgc_2022_cluster_map_v2",
    gene_specificity_matrix_source=SpecificityMatrixSource(
        GTEx_V10_MEDIAN_TISSUE_EXPRESSION_RNA_SEQ_PREP_FOR_MAGMA,
        pipe=DropColPipe(["Average"]),
        gene_col="Gene",
    ),
    extracted_magma_gene_results_source=GeneInfoSource(
        SCH_PGC_2022_STANDARD_ANALYSIS.magma_tasks.inner.gene_analysis_extracted_result_task,
        pipe=IdentityPipe(),
        gene_col=MAGMA_GENE_COL,
    ),
    gene_thesaurus_source=GeneInfoSource(
        GENE_THESAURUS, pipe=IdentityPipe(), gene_col="Gene stable ID"
    ),
)
