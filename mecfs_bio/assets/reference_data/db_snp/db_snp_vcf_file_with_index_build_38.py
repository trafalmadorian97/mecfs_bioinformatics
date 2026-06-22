"""
Task that downloads a vcf file containing the genome build 38 dbSNP reference database into a directory,
then indexes it using bcftools.

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.download_files_into_directory_task import (
    DownloadEntry,
    DownloadFilesIntoDirectoryTask,
)

DB_SNP_VCF_FILE_WITH_INDEX_BUILD_38_DIR = DownloadFilesIntoDirectoryTask(
    meta=ReferenceDataDirectoryMeta(
        id="db_snp_build_38_with_index",
        group="db_snp_reference_data",
        sub_group="build_38",
        sub_folder=PurePath("vcf"),
    ),
    entries=[
        DownloadEntry.create_with_bcftools_index(
            url="https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz",
            md5hash=None,
            filename="GCF_000001405.40.gz",
        )
    ],
)
