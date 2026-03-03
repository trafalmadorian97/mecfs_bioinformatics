"""
Summary statistics from Ankylosing spondylitis in Finngen


See: https://opengwas.io/datasets/finn-b-M13_ANKYLOSPON#

opengwas only provides a temporary link, so I mirrored on dropbox.
"""
from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
import polars as pl

FINNGEN_ANKYLOSING_SPONDYLITIS_DATA_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("finngen_spond_eur_raw"),
        trait="ankylosing_spondylitis",
        project="finngne",
        sub_dir="raw",
        project_path=PurePath("finn-b-M13_ANKYLOSPON.vcf.gz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t",
                                                               has_header=False,
                                                               column_names=["CHROM",
                                                                             "POS",
                                                                             "ID",
                                                                             "REF",
                                                                             "ALT",
                                                                             "QUAL",
                                                                             "FILTER",
                                                                             "INFO",
                                                                             "FORMAT",
                                                                             "finn-b-M13_ANKYLOSPON"
                                                                             ],
                                                               comment_char="#",
                                                               schema_overrides={
                                                                   "CHROM":pl.String()
                                                               }

                                                               )),
    ),
    url="https://www.dropbox.com/scl/fi/keanskgpwrf26ngdly3dm/finn-b-M13_ANKYLOSPON.vcf.gz?rlkey=gwzualqvpiyktd965yhpdpztn&dl=1",
    md5_hash=None,
)
