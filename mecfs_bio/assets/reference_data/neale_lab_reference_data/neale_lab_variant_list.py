from pathlib import PurePath

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

NEALE_LAB_VARIANTS_REFERENCE = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="neale_lab",
        sub_group="metadata",
        sub_folder=PurePath("variants"),
        extension=".bgz",
        filename="variants.tsv",
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat("\t", schema_overrides={"chr": pl.String()})
        ),
        id=AssetId("neale_lab_variants_reference"),
    ),
    url="https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz",
    md5_hash="4fc9936ba4b4fd446cb4325dd63b6e72",
)
