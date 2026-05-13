from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

NEALE_LAB_PHENOTYPE_SUMMARY_FILE = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="neale_lab",
        sub_group="metadata",
        sub_folder=PurePath("phenotypes"),
        extension=".bgz",
        filename="phenotypes.both_sexes.tsv",
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat("\t")),
        id=AssetId("neale_lab_phenotypes_reference"),
    ),
    url=" https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/phenotypes.both_sexes.tsv.bgz",
    md5_hash="cba910ee6f93eaed9d318edcd3f1ce18",
)
