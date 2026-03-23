"""
Raw Summary Statistics from the [Million Veterans Program](https://www.science.org/doi/abs/10.1126/science.adj1182) GWAS of myocardial infraction


See here for the ftp download page: https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475932/harmonised/

See the GWAS catalog page here: https://www.ebi.ac.uk/gwas/studies/GCST90475932

The causal chain LDL -> Atherosclerosis -> Myocardial Infarction is well-established, and thus can serve as
a positive control for the detection of causal pathways.  That is, any causal pathway detection method should
be able to detect this pathway.


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERAN_MI_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("million_veterans_mi_eur_raw"),
        trait="myocardial_infaction",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90475932.tsv.gzip"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475932/harmonised/GCST90475932.h.tsv.gz",
    md5_hash="ddf8704acb4732834013fd907dc9fa60",
)
