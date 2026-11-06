"""
Task to download van Rheenen et al. 2021's amyotrophic lateral sclerosis GWAS, European
ancestry (27,205 cases, 110,881 controls), as deposited on the GWAS Catalog (GCST90027164).

This is the submitted GRCh37 file, not the Catalog's harmonised one: only the submitted file has
the per-variant N_effective column. Alleles are lowercase.
See https://www.ebi.ac.uk/gwas/studies/GCST90027164

Citation:
van Rheenen, Wouter, et al. "Common and rare variant association analyses in amyotrophic
lateral sclerosis identify 15 risk loci with distinct genetic architectures and
neuron-specific biology." Nature Genetics 53.12 (2021): 1636-1648.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

VAN_RHEENEN_2021_ALS_EUR_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("van_rheenen_2021_als_eur_raw"),
        trait="als",
        project="van_rheenen_et_al_2021",
        sub_dir="raw",
        project_path=PurePath("GCST90027164_buildGRCh37.tsv.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(separator="\t", null_values=["NaN", "NA"]),
        ),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz",
    md5_hash="705fa6657d75163d313138798ce9a579",
)
