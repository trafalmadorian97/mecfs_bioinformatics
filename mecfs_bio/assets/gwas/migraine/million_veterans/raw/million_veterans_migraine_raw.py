"""
This Task downloads the raw summary statistics data for the million-veterans GWAS of migraines from
GWAS catalog.

For GWAS catalog study info see: https://www.ebi.ac.uk/gwas/studies/GCST90475837


Citation: Verma, Anurag, et al. "Diversity and scale:
Genetic architecture of 2068 traits in the VA Million Veteran Program." Science 385.6706 (2024): eadj1182.


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERAN_MIGRAINE_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("million_veterans_migraine_eur_raw"),
        trait="migraine",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90475837.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475837/harmonised/GCST90475837.h.tsv.gz",
    md5_hash="b5b819a748fae3555838fd55328c41ae",
)
