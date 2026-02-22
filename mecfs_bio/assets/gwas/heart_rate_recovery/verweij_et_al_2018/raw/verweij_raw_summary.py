"""
GWAS summary statistics from Verweij, Niek, Yordi J. van de Vegte, and Pim van der Harst. "Genetic study links components of the autonomous nervous system to heart-rate profile during exercise." Nature communications 9.1 (2018): 898.

Download from GWAS Catalog ftp site.

See also: https://www.ebi.ac.uk/gwas/studies/GCST005850

"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

VERWEIJI_ET_AL_RAW_HARMONIZED_BUILD_37 = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("verweij_et_al_2018_raw_harmonized"),
        trait="heart_rate_recovery",
        project="verweij_et_al_2018",
        sub_dir="raw",
        project_path=PurePath("29497042-GCST005850-EFO_0009185-build37.f.tsv.gz"),
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
    ),
    url="http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST005001-GCST006000/GCST005850/harmonised/29497042-GCST005850-EFO_0009185-build37.f.tsv.gz",
    md5_hash="076ec2bf05476ce518fb628ad0bec8f2",
)
