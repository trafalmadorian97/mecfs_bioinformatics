from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.magma.decode_me_gwas_1_build_37_magma_snp_locs import (
    DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_LOCS,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_37_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED,
)
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask

DECODE_ME_GWAS_1_MAGMA_ENTREZ_ANNOTATIONS = MagmaAnnotateTask.create(
    asset_id="decode_me_gwas_1_build_37_magma_annotations",
    magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
    snp_loc_file_task=DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_LOCS,
    gene_loc_file_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED,
)
