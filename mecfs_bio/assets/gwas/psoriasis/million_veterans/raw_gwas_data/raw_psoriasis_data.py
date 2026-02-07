"""
Raw Summary Statistics from the [Million Veterans Program](https://www.science.org/doi/abs/10.1126/science.adj1182) GWAS of Psoriasis.


See here for the ftp download page: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476186/

See the GWAS catalog page here: https://www.ebi.ac.uk/gwas/studies/GCST90476186

Note that this is a gwas of Psoriasis and related disorder (PheCode 696).


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERAN_PSORIASIS_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("million_veterans_psoriasis_eur_raw"),
        trait="psoriasis",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90476186.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90476001-GCST90477000/GCST90476186/GCST90476186.tsv.gz",
    md5_hash="72c556837df9e1d3df092a484d3a4449",
)
