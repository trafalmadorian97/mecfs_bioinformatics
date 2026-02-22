"""
Raw Summary Statistics from the [Million Veterans Program](https://www.science.org/doi/abs/10.1126/science.adj1182) GWAS of LDL cholesterol.


See here for the ftp download page:https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475416/

See the GWAS catalog page here: https://www.ebi.ac.uk/gwas/studies/GCST90475416


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERAN_LDL_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("million_veterans_ldl_eur_raw"),
        trait="ldl",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90475416.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475416/GCST90475416.tsv.gz",
    md5_hash="ecc667e9619db80b1752299bd81eaac5",
)
