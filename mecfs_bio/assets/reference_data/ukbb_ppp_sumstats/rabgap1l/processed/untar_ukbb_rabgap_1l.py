from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.raw.ukbb_ppp_rabgap1l import (
    UKBB_PPP_RABGAP1L,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

UKBB_PPP_RABGAP1L_UNTAR = ExtractTarGzipTask.create_from_reference_file_task(
    asset_id="ukbb_ppp_rabgap_1l_untar", source_task=UKBB_PPP_RABGAP1L, read_mode="r"
)
