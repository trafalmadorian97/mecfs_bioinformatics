"""
Task to download GWAS summary statistics from Johnston et al's multisite pain GWAS

see:
https://www.ebi.ac.uk/gwas/studies/GCST008512


Citation:
Johnston, Keira J A et al. “Genome-wide association study of multisite chronic pain in UK Biobank.” PLoS genetics vol. 15,6 e1008164. 13 Jun. 2019, doi:10.1371/journal.pgen.1008164
"""

from pathlib import PurePath

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

JOHNSTON_ET_AL_PAIN_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("johnston_et_al_pain_raw"),
        trait="multisite_pain",
        project="johnston_et_al",
        sub_dir="raw",
        project_path=PurePath("chronic_pain-bgen.stats.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(
                separator="\t",
                null_values=["NaN", "NA"],
                schema_overrides={"GENPOS": pl.Float32()},
            ),
        ),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST008001-GCST009000/GCST008512/chronic_pain-bgen.stats.gz",
    md5_hash="cc175d6a6e5ee80005f3bdb78b0b7694",
)
