"""
GWAS summary statistics from the paper: Lee, James J., et al. "Gene discovery and polygenic prediction from a genome-wide association study of educational attainment in 1.1 million individuals." Nature genetics 50.8 (2018): 1112-1121.

Harmonized to fit the GWAS catalog format.

GWAS catalog page is here: https://www.ebi.ac.uk/gwas/publications/30038396.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

LEE_ET_AL_2018_CATALOG_HARMONIZED_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("lee_at_al_2018_catalog_harmonized"),
        trait="educational_attainment",
        project="lee_el_al_2018",
        sub_dir="raw",
        project_path=PurePath(
            "harmonized/30038396-GCST006572-EFO_0008354-build37.f.tsv.gz"
        ),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST006001-GCST007000/GCST006572/harmonised/30038396-GCST006572-EFO_0008354-build37.f.tsv.gz",
    md5_hash=None,
)
