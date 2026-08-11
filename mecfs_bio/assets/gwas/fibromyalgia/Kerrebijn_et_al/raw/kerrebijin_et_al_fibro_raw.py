"""
Download raw European summary statistics from kerrebijin et al's GWAS of fibromyalgia

Citation:
Kerrebijn, Isabel, et al. "The genetic architecture of fibromyalgia across 2.5 million individuals." medRxiv (2025).


See GWAS Catalog entry: https://www.ebi.ac.uk/gwas/studies/GCST90838603

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

KERREBIJN_ET_AL_FIBRO_EUR_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("kerrebijn_et_al_eur_fibro_raw"),
        trait="fibromyalgia",
        project="kerrebijn_et_al",
        sub_dir="raw",
        project_path=PurePath("GCST90838603.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90838001-GCST90839000/GCST90838603/harmonised/GCST90838603.h.tsv.gz",
    md5_hash="3e2ea8a5e951f1810f9c23ea01573e44",
)
