from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
# see: https://cloufield.github.io/gwaslab/AssignrsID/
"""
DB_SNP_VCF_FILE_BUILD_37 = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="db_snp_reference_data",
        sub_group="build_37",
        sub_folder=PurePath("vcf"),
        asset_id="db_snp_vcf_reference_file_build_37",
        filename="GCF_000001405.25",
        extension=".gz",
    ),
    url="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz",
    md5_hash="35db22bcd166f904e4775dbbc29f5965",
)
