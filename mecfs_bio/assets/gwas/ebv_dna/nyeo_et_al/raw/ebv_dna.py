"""
Task to download raw summary statistics from Nyeo et al's GWAS of EBV DNA.

Citation: Nyeo, Sherry S., et al. "Population-scale sequencing resolves determinants of persistent EBV DNA." Nature (2026): 1-9.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

NYEO_EBV_DNA_SUMSTATS = DownloadFileTask(
    url="https://my.locuszoom.org/gwas/409414/data/",
    meta=GWASSummaryDataFileMeta(
        id="nyeo_ebv_dna_sumstats",
        trait="ebv_dna",
        project="nyeo_et_al",
        sub_dir="raw",
        project_path=PurePath("ebv_dna_sumstats.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    md5_hash="5bc20da4b176fe1c605f35fd69577783",
)
