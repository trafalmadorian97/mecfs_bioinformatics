from mecfs_bio.assets.reference_data.lava_ld_reference.ukbb.raw.ukbb_lava_ld_ref import (
    LAVA_UKBB_LD_REF,
)
from mecfs_bio.build_system.task.extract_all_zips_in_directory_task import (
    ExtractAllZipsInDir,
)

LAVA_UKBB_LD_REF_EXTRACTED = ExtractAllZipsInDir.create(
    source_directory_task=LAVA_UKBB_LD_REF, asset_id="lava_ukbb_reference_ld_extracted"
)
