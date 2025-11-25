"""
Pull raw gwas summary data from:
"Meta-analysis of genome-wide association studies for body fat distribution in 694,649 individuals of European ancestry." Pulit, SL et al. bioRxiv, 2018. https://www.biorxiv.org/content/early/2018/04/18/304030
Files are hosted at: https://zenodo.org/records/1251813
"""
from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameTextFormat, DataFrameReadSpec
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PULLIT_ET_AL_2018_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        short_id=AssetId("pullit_et_al_2018_gwas_summary"),
        trait="body_fat_distribution",
        project="pullt_et_al_2018",
        sub_dir="raw",
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator=" ")),
        project_path=PurePath("whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz"),
    ),
    url="https://zenodo.org/records/1251813/files/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz?download=1",
    md5_hash="500ee10b8406dc80394ad8cc968a4eb1",
)
