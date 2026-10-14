from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.magma.liu_et_al_2023_eur_37_magma_snp_locs import (
    LIU_ET_AL_IBD_EUR_BUILD_37_MAGMA_SNP_LOCS,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask

LIU_ET_AL_2023_EUR_IBD_37_ENSEMBL_ANNOTATIONS = MagmaAnnotateTask.create(
    asset_id="liu_et_al_2023_eur_37_magma_ensemble_annotations",
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    snp_loc_file_task=LIU_ET_AL_IBD_EUR_BUILD_37_MAGMA_SNP_LOCS,
    gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
