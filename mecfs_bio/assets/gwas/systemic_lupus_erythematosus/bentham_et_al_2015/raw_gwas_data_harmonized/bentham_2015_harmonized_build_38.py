from pathlib import PurePath

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

"""
See: https://www.ebi.ac.uk/gwas/publications/26502338
for the gwas catalog entry for this publication
"""
BENTHAM_2015_HARMONIZED_BUILD_38 = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("bentham_2015_harmonized_build_38"),
        trait="systemic_lupus_erythematosus",
        project="bentham_et_al_2015",
        project_path=PurePath("26502338-GCST003156-EFO_0002690.h.tsv.gz"),
        sub_dir="raw_harmonized",
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(
                separator="\t",
                null_values=["NA"],
                schema_overrides={"hm_chrom": pl.String()},
            )
        ),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/harmonised/26502338-GCST003156-EFO_0002690.h.tsv.gz",
    md5_hash="bae52c67bcea4b8e4308eb6f0c45f924",
)
