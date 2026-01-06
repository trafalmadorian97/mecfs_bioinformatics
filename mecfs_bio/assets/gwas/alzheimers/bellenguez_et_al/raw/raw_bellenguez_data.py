"""
Task to download GWAS summary statistics from Bellenguez et al's Alzhiemers GWAS

see:
https://www.ebi.ac.uk/gwas/studies/GCST90027158


Citation:
Bellenguez, Céline, et al. "New insights into the genetic etiology of Alzheimer’s disease and related dementias.Nature genetics (2022): 412-436.
"""
from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

BELLENGUEZ_ET_AL_ALZHIEMERS_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        short_id=AssetId("bellenguez_et_al_alzhiemers_raw"),
        trait="alzheimers",
        project="bellenguez_et_al",
        sub_dir="raw",
        project_path=PurePath("35379992-GCST90027158-MONDO_0004975.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027158/harmonised/35379992-GCST90027158-MONDO_0004975.h.tsv.gz",
    md5_hash="e2c8be73b5aa7698c5e8878ae607fe85",
)
