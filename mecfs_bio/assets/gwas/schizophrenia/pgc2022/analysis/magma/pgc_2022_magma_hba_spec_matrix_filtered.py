from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR,
)
from mecfs_bio.assets.reference_data.magma_specificity_matrices.raw.magma_specificity_matrix_from_hbca_rna_duncan import (
    MAGMA_ENTREZ_SPECIFICITY_MATRIX_HBCA_RNA_DUNCAN,
)
from mecfs_bio.build_system.task.magma.magma_subset_specificity_matrix_using_top_labels import (
    MagmaSubsetSpecificityMatrixWithTopLabels,
)

PGC2022_HBA_MAGMA_SPEC_MATRIX_FILTERED = (
    MagmaSubsetSpecificityMatrixWithTopLabels.create(
        asset_id="pgc2022_magma_hba_spec_matrix_filtered",
        specificity_matrix_task=MAGMA_ENTREZ_SPECIFICITY_MATRIX_HBCA_RNA_DUNCAN,
        magma_gene_covar_analysis_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR,
    )
)
