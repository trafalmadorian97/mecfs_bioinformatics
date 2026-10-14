from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_raw import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_REF,
)
from mecfs_bio.build_system.task.extract_all_from_zip_task import ExtractAllFromZipTask

MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED = (
    ExtractAllFromZipTask.create_from_zipped_reference_file(
        source_task=MAGMA_EUR_BUILD_37_1K_GENOMES_REF,
        asset_id="magma_eur_1k_genomes_build_37_ld_ref_extracted",
    )
)
