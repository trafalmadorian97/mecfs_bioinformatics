"""
Task that downloads a vcf file describing a reference database of SNPs into a directory,
then indexes it using bcftools

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.download_files_into_directory_task import (
    DownloadEntry,
    DownloadFilesIntoDirectoryTask,
)

DB_SNP_VCF_FILE_WITH_INDEX_BUILD_37_DIR = DownloadFilesIntoDirectoryTask(
    meta=ReferenceDataDirectoryMeta(
        asset_id="db_snp_build_37_with_index",
        group="db_snp_reference_data",
        sub_group="build_37",
        sub_folder=PurePath("vcf"),
    ),
    entries=[
        DownloadEntry.create_with_bcftools_index(
            url="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz",
            md5hash="35db22bcd166f904e4775dbbc29f5965",
            filename="GCF_000001405.25.gz",
        )
    ],
)
