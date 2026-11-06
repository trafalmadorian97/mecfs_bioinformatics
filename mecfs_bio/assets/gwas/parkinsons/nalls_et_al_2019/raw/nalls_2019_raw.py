"""
Task to download the GWAS Catalog harmonised summary statistics of Nalls et al. 2019's
Parkinson's disease meta-GWAS without 23andMe (GCST009325): 15,056 cases and 18,618 UK Biobank
proxy cases (a parent or sibling with Parkinson's) against 449,056 controls, European ancestry.

The harmonised file is GRCh38 (lifted over by the Catalog) and carries per-variant N_cases and
N_controls. See https://www.ebi.ac.uk/gwas/studies/GCST009325

Citation:
Nalls, Mike A., et al. "Identification of novel risk loci, causal insights, and heritable risk
for Parkinson's disease: a meta-analysis of genome-wide association studies." The Lancet
Neurology 18.12 (2019): 1091-1102.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

NALLS_2019_PARKINSONS_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("nalls_2019_parkinsons_raw"),
        trait="parkinsons",
        project="nalls_et_al_2019",
        sub_dir="raw",
        project_path=PurePath("GCST009325.h.tsv.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(separator="\t", null_values=["NaN", "NA"]),
        ),
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz",
    md5_hash="e7f2654a34931bb9996fb0927cece282",
    # The EBI FTP server limits speed per connection (~200 KB/s observed).
    connections=16,
)
