"""
This Task downloads the raw summary statistics data for the UK Biobank GWAS of migraines from
GWAS catalog.

For GWAS catalog study info see: https://www.ebi.ac.uk/gwas/studies/GCST90473326


Citation: UK Biobank Whole-Genome Sequencing Consortium.
“Whole-genome sequencing of 490,640 UK Biobank participants.”
Nature vol. 645,8081 (2025): 692-701. doi:10.1038/s41586-025-09272-9


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

UK_BIOBANK_2025_MIGRAINE_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("uk_biobank_2025_migraine_eur_raw"),
        trait="migraine",
        project="uk_biobank_2025",
        sub_dir="raw",
        project_path=PurePath("GCST90473326.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90473001-GCST90474000/GCST90473326/harmonised/GCST90473326.h.tsv.gz",
    md5_hash="3ca307936baa3564cbfb4f95e04bbc1d",
)
