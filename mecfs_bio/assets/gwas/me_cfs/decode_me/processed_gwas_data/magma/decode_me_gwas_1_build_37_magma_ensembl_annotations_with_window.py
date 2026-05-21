from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_snp_locs import (
    DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_LOCS,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask

DECODE_ME_GWAS_1_MAGMA_ENSEMBL_ANNOTATIONS_WITH_WINDOW = MagmaAnnotateTask.create(
    asset_id="decode_me_gwas_1_build_37_magma_ensembl_annotations_with_window",
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    snp_loc_file_task=DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_LOCS,
    gene_loc_file_task=MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    window=(35, 10),  # This choice comes from the Duncan paper
)
