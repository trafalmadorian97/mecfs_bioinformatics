"""
Task to download raw summary statistics of Han et al's asthma GWAS as hosted on the GWAS catalog

Paper: Han, Yi, et al. "Genome-wide analysis highlights contribution of immune system pathways to the genetic architecture of asthma." Nature communications 11.1 (2020): 1776.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

HAN_ET_AL_ASTHMA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("han_et_al_asthma_raw"),
        trait="asthma",
        project="han_et_al",
        sub_dir="raw",
        project_path=PurePath("HanY_prePMID_asthma_UKBB.txt.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST010001-GCST011000/GCST010042/HanY_prePMID_asthma_UKBB.txt.gz",
    md5_hash="0d34c6cb3581e5d8af73797bfa256f68",
)
