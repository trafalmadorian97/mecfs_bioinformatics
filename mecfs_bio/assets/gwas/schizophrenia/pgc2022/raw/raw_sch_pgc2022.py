# https://figshare.com/ndownloader/files/34517828
from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PGC_2022_SCH_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        short_id=AssetId("pgc_2022_sch_raw"),
        trait="schizophrenia",
        project="pgc_2022",
        sub_dir="raw",
        project_path=PurePath("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(separator="\t", comment_char="#")
        ),
    ),
    url="https://figshare.com/ndownloader/files/34517828",
    md5_hash="6ebe2376f5cda972d37efa0f214c4df0",
)
