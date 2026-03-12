from mecfs_bio.assets.reference_data.mixer.raw.mixer_g1000_plink_eur_raw import MIXER_RAW_G1000_PLINK_DATA
from mecfs_bio.build_system.task.extract_all_from_zip_task import ExtractAllFromZipTask, UseCommandLineTool

MIXER_G1000_PLINK_EXTRACTED=ExtractAllFromZipTask.create_from_zipped_reference_file(
    source_task=MIXER_RAW_G1000_PLINK_DATA,
    asset_id="mixer_g1000_plink_eur_extracted",
    # use_command_line=UseCommandLineTool()
)