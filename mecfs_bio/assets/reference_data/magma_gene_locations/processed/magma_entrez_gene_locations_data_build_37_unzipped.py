from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_entrez_gene_location_reference_data_build_37 import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import (
    ExtractFromZipTask,
)

MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_EXTRACTED = (
    ExtractFromZipTask.create_from_zipped_reference_file(
        source_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
        asset_id="magma_gene_location_reference_data_build_37_extracted",
        file_to_extract="NCBI37.3.gene.loc",
    )
)
