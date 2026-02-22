from pathlib import PurePath

from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.reference.schemas.hg19_sn151_schema import (
    HG19_SNP151_SCHEMA,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
See
https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/
"""
GENOME_ANNOTATION_DATABASE_BUILD_37 = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="genome_annotations",
        sub_group="build_37",
        sub_folder=PurePath("raw"),
        id="genome_annotation_database_build_37",
        extension=".gz",
        filename="snp151.txt",
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(
                separator="\t",
                column_names=HG19_SNP151_SCHEMA,
                has_header=False,
            )
        ),
    ),
    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151.txt.gz",
    md5_hash=None,
)
