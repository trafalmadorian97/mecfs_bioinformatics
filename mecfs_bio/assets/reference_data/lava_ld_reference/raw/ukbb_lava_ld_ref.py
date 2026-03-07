from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import ReferenceDataDirectoryMeta
from mecfs_bio.build_system.task.download_files_into_directory_task import DownloadFilesIntoDirectoryTask, DownloadEntry

#
#
#+ wget https://vu.data.surfsara.nl/index.php/s/7NBVIvtPRdu7Qhz/download -O lava-ukb-v1.1_chr1-2.zip
#+ wget https://vu.data.surfsara.nl/index.php/s/fy6ITboMojHrQXr/download -O lava-ukb-v1.1_chr3-4.zip
#+ wget https://vu.data.surfsara.nl/index.php/s/mRz31q0lq7KMcuI/download -O lava-ukb-v1.1_chr5-6.zip
#+ wget https://vu.data.surfsara.nl/index.php/s/89RXxPN2BlOqxLs/download -O lava-ukb-v1.1_chr7-9.zip
#+ wget https://vu.data.surfsara.nl/index.php/s/3dU1L2Hap43xuCs/download -O lava-ukb-v1.1_chr10-12.zip
# wget https://vu.data.surfsara.nl/index.php/s/DQqJ2Sqr49RP4xe/download -O lava-ukb-v1.1_chr13-16.zip
# wget https://vu.data.surfsara.nl/index.php/s/U8eg5XfTr8qPzPp/download -O lava-ukb-v1.1_chr17-23.zip
LAVA_LD_REF=DownloadFilesIntoDirectoryTask(
    meta=ReferenceDataDirectoryMeta(
        group="lavareference_ld",
        sub_group="ukbb",
        sub_folder=PurePath("raw"),
        id=AssetId("lava_ukbb_reference_ld")
    ),
    entries=[
        DownloadEntry(
           url= "https://vu.data.surfsara.nl/index.php/s/7NBVIvtPRdu7Qhz/download",
            filename="lava-ukb-v1.1_chr1-2.zip",
            md5hash=None,
            post_download_command=None
        ),

        DownloadEntry(
            url="https://vu.data.surfsara.nl/index.php/s/fy6ITboMojHrQXr/download",
            filename="lava-ukb-v1.1_chr3-4.zip",
            md5hash=None,
            post_download_command=None
        ),

        DownloadEntry(
            url="https://vu.data.surfsara.nl/index.php/s/mRz31q0lq7KMcuI/download",
            filename="lava-ukb-v1.1_chr5-6.zip",
            md5hash=None,
            post_download_command=None
        ),

        DownloadEntry(
            url="https://vu.data.surfsara.nl/index.php/s/89RXxPN2BlOqxLs/download",
            filename="lava-ukb-v1.1_chr7-9.zip",
            md5hash=None,
            post_download_command=None
        ),

        DownloadEntry(
            url="https://vu.data.surfsara.nl/index.php/s/3dU1L2Hap43xuCs/download",
            filename="lava-ukb-v1.1_chr10-12.zip",
            md5hash=None,
            post_download_command=None
        ),

        DownloadEntry(
            url="https://vu.data.surfsara.nl/index.php/s/DQqJ2Sqr49RP4xe/download",
            filename="lava-ukb-v1.1_chr13-16.zip",
            md5hash=None,
            post_download_command=None
        ),

        DownloadEntry(
            url="https://vu.data.surfsara.nl/index.php/s/U8eg5XfTr8qPzPp/download",
            filename="lava-ukb-v1.1_chr17-23.zip",
            md5hash=None,
            post_download_command=None
        ),

    ]
)