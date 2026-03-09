"""
See: https://www.ebi.ac.uk/gwas/publications/26502338
for the gwas catalog entry for this publication

Summary statistics linked from here: https://www.ebi.ac.uk/gwas/studies/GCST003156
According to the Gwas catalog metadata, the sample size is 14267

Citation: Bentham, James, et al. "Genetic association analyses implicate aberrant regulation of innate and adaptive
 immunity genes in the pathogenesis of systemic lupus erythematosus." Nature genetics 47.12 (2015): 1457-1464.

"""

import polars as pl

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

BENTHAM_2015_RAW_BUILD_37 = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("bentham_2015_raw_build_37"),
        trait="systemic_lupus_erythematosus",
        project="bentham_et_al_2015",
        sub_dir="raw",
        project_path=None,
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(
                separator="\t",
                null_values=["NA"],
                schema_overrides={"chrom": pl.String()},
            )
        ),
        extension=".tsv.gz",
    ),
    url="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST003001-GCST004000/GCST003156/bentham_2015_26502338_sle_efo0002690_1_gwas.sumstats.tsv.gz",
    md5_hash=None,
)
