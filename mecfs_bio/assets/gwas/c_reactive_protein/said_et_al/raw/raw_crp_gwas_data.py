from pathlib import PurePath

import polars

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

SAID_CRP_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("said_crp_eur_raw"),
        trait="crp",
        project="said_et_al_2022",
        sub_dir="raw",
        project_path=PurePath("35459240-GCST90029070-EFO_0004458.h.tsv.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(
                separator="\t",
                null_values=["NA"],
                schema_overrides={"hm_chrom": polars.String()},
            )
        ),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90029001-GCST90030000/GCST90029070/harmonised/35459240-GCST90029070-EFO_0004458.h.tsv.gz",
    md5_hash="8156fe76b6149c33cb31232cd949bc1e",
)
