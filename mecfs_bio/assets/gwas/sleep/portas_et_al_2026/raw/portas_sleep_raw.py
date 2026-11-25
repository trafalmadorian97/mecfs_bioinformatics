"""
Tasks to download GWAS Catalog harmonised summary statistics for the device-measured sleep traits from:

Portas, Laura, et al. "Genetic architecture of sleep in a genome wide association study of device measured sleep traits." Nature Communications (2026).

Publication page: https://www.ebi.ac.uk/gwas/publications/41922918

Sleep traits were derived from UK Biobank wrist-worn accelerometer recordings by a deep learning model trained on
polysomnography data.  Association analyses used REGENIE on 80,013 European ancestry participants (sexes combined).
The authors excluded variants with MAF < 0.01 or INFO < 0.3.

Only the sexes-combined studies have full summary statistics in the GWAS Catalog; the sex-stratified studies do not.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

PORTAS_SLEEP_DURATION_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("portas_et_al_2026_sleep_duration_raw"),
        trait="sleep_duration",
        project="portas_et_al_2026",
        sub_dir="raw",
        project_path=PurePath("GCST90824104.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/GCST90824104/harmonised/GCST90824104.h.tsv.gz",
    md5_hash="8aad7906354bf7c50034d864e5d44a73",
)

PORTAS_SLEEP_EFFICIENCY_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("portas_et_al_2026_sleep_efficiency_raw"),
        trait="sleep_efficiency",
        project="portas_et_al_2026",
        sub_dir="raw",
        project_path=PurePath("GCST90824105.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/GCST90824105/harmonised/GCST90824105.h.tsv.gz",
    md5_hash="ea714c140248e110319f37f35c886ceb",
)

PORTAS_REM_SLEEP_DURATION_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("portas_et_al_2026_rem_sleep_duration_raw"),
        trait="rem_sleep_duration",
        project="portas_et_al_2026",
        sub_dir="raw",
        project_path=PurePath("GCST90824106.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/GCST90824106/harmonised/GCST90824106.h.tsv.gz",
    md5_hash="21e7cf759ff1656cf54c68b1d60231fc",
)

PORTAS_NREM_SLEEP_DURATION_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("portas_et_al_2026_nrem_sleep_duration_raw"),
        trait="nrem_sleep_duration",
        project="portas_et_al_2026",
        sub_dir="raw",
        project_path=PurePath("GCST90824107.h.tsv.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90824001-GCST90825000/GCST90824107/harmonised/GCST90824107.h.tsv.gz",
    md5_hash="be1437276bce5c9678ca93f614dc68cf",
)
