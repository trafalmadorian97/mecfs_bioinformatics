"""
This Task downloads the raw summary statistics data for the FinnGen release 13 GWAS of migraine with aura
(endpoint G6_MIGRAINE_WITH_AURA: 12,743 cases, 366,556 controls).

Coordinates are on GRCh38.

For FinnGen endpoint info see: https://risteys.finngen.fi/endpoints/G6_MIGRAINE_WITH_AURA


Citation: Kurki, Mitja I., et al. "FinnGen provides genetic insights from a well-phenotyped
isolated population." Nature 613.7944 (2023): 508-518. doi:10.1038/s41586-022-05473-8


"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

FINNGEN_R13_MIGRAINE_WITH_AURA_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("finngen_r13_migraine_with_aura_raw"),
        trait="migraine",
        project="finngen_r13_migraine_with_aura",
        sub_dir="raw",
        project_path=PurePath("finngen_R13_G6_MIGRAINE_WITH_AURA.gz"),
        read_spec=DataFrameReadSpec(
            format=DataFrameTextFormat(separator="\t", null_values=["NA"])
        ),
    ),
    url="https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_G6_MIGRAINE_WITH_AURA.gz",
    md5_hash="abb7380162e534d43b436dc59fd9dbe8",
)
