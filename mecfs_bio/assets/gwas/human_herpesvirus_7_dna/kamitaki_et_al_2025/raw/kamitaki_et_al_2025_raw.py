"""
Script to download summary statistics from kamitaki et al'2 GWAS of Human Herpesvirus 7 DNA

Citation:
Kamitaki, Nolan, et al. "Genes and environment profoundly affect the human virome." bioRxiv (2025): 2025-09.

link:https://www.biorxiv.org/content/10.1101/2025.09.08.674901
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

KAMITAKI_ET_AL_HHV7_DNA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id="kamitaki_et_al_2025_hhb7_raw",
        trait="hhb7_dna",
        project="kamataki_et_al",
        sub_dir="raw",
        project_path=PurePath("HHV7_invnorm.bgen.MAF.stats.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://data.broadinstitute.org/lohlab/virome_summary_statistics/blood/HHV7_invnorm.bgen.MAF.stats.gz",
    md5_hash="e4a7a0b36f23413d9de7da4bb2b66c18",
)
