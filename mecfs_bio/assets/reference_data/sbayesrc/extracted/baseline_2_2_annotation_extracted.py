"""
SBayesRC baseline-LD model 2.2 functional annotation (extracted file).

Extracts the single annot_baseline2.2.txt member directly to a FileAsset, suitable
for passing as SBayesRCTask's annotation_file_task.
"""

from mecfs_bio.assets.reference_data.sbayesrc.raw.baseline_2_2_annotation_raw import (
    SBAYESRC_BASELINE_2_2_ANNOTATION_RAW,
)
from mecfs_bio.build_system.task.extraction_one_file_from_zip_task import (
    ExtractFromZipTask,
)

SBAYESRC_BASELINE_2_2_ANNOTATION_EXTRACTED = (
    ExtractFromZipTask.create_from_zipped_reference_file(
        source_task=SBAYESRC_BASELINE_2_2_ANNOTATION_RAW,
        asset_id="sbayesrc_baseline_2_2_annotation_extracted",
        file_to_extract="annot_baseline2.2.txt",
        extension="",
    )
)
