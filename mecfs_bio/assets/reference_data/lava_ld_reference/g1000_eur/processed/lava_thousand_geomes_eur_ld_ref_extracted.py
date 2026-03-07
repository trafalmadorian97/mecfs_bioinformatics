from mecfs_bio.assets.reference_data.lava_ld_reference.g1000_eur.raw.thousand_genomes_eur_ld_ref import \
    LAVA_G1000_EUR_REF_LD
from mecfs_bio.build_system.task.extract_all_from_zip_task import ExtractAllFromZipTask

LAVA_G100_EUR_LD_REF_EXTRACTED = ExtractAllFromZipTask.create_from_zipped_reference_file(
LAVA_G1000_EUR_REF_LD,
    asset_id="lava_g1000_eur_ld_ref_extracted"
)