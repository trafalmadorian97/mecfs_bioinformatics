"""

Raw Summary Statistics from the [Million Veterans Program](https://www.science.org/doi/abs/10.1126/science.adj1182) GWAS of red-blood cell volume amomng Europeans.

See:

https://www.ebi.ac.uk/gwas/studies/GCST90475466
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        short_id=AssetId("million_veterans_red_blood_cell_volume_eur_raw"),
        trait="red_blood_cell_volume",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90475466.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475466/GCST90475466.tsv.gz",
    md5_hash="62f0eda1fdac54a06b8259f137a43ddd",
)

MILLION_VETERAN_RED_BLOOD_VOLUME_EUR_DATA_RAW_HARMONIZED = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        short_id=AssetId("million_veterans_red_blood_cell_volume_eur_raw_harmonized"),
        trait="red_blood_cell_volume",
        project="million_veterans",
        sub_dir="raw",
        project_path=PurePath("GCST90475466.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90475001-GCST90476000/GCST90475466/harmonised/GCST90475466.h.tsv.gz",
    md5_hash="3e96c443607751e95419359dde49d3b6",
)
