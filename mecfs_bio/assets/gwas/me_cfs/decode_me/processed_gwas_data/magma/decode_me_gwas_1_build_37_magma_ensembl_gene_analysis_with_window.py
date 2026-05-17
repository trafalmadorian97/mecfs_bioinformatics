from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_ensembl_annotations_with_window import (
    DECODE_ME_GWAS_1_MAGMA_ENSEMBL_ANNOTATIONS_WITH_WINDOW,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_snp_p_values import (
    DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_P_VALUES,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)

DECODE_ME_GWAS_1_MAGMA_ENSEMBL_GENE_ANALYSIS_WITH_WINDOW = MagmaGeneAnalysisTask.create(
    asset_id="decode_me_gwas_1_build_37_magma_ensemble_gene_analysis_with_window",
    magma_annotation_task=DECODE_ME_GWAS_1_MAGMA_ENSEMBL_ANNOTATIONS_WITH_WINDOW,
    magma_p_value_task=DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_P_VALUES,
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
    ld_ref_file_stem="g1000_eur",
    sample_size=275488,  # source: add up cases and controls column in raw gwas data file
)
