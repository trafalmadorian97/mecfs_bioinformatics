"""
SBayesRC UKB EUR HapMap3 LD reference (extracted).

Extracts the ukbEUR_HM3 directory from the downloaded archive; the resulting
DirectoryAsset is what SBayesRCTask and the polypwas tasks pass as their LD
reference.
"""

from mecfs_bio.assets.reference_data.sbayesrc.raw.ukb_eur_hm3_ld_raw import (
    SBAYESRC_UKB_EUR_HM3_LD_RAW,
)
from mecfs_bio.build_system.task.extract_all_from_zip_task import ExtractAllFromZipTask

SBAYESRC_UKB_EUR_HM3_LD_EXTRACTED = (
    ExtractAllFromZipTask.create_from_zipped_reference_file(
        source_task=SBAYESRC_UKB_EUR_HM3_LD_RAW,
        asset_id="sbayesrc_ukb_eur_hm3_ld_extracted",
        sub_folder_name_inside_zip="ukbEUR_HM3",
    )
)
