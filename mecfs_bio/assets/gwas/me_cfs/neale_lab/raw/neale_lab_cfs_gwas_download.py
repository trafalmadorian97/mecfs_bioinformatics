"""
This task downloads the GWAS summary statistics from UK Biobank GWAS of self-reported CFS performed by the
Neale Lab


The association model used by the Neale lab is described on the blog post here:
http://www.nealelab.is/blog/2017/9/11/details-and-considerations-of-the-uk-biobank-gwas

Notably, for efficiency the authors used linear regression instead of logistic regression even for case-control stduies.

A Google sheet that serves as an index of Neale Lab GWAS can be found here:
https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit?gid=178908679#gid=178908679`
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

NEALE_LAB_CFS_RAW = DownloadFileTask(
    meta=GWASSummaryDataFileMeta(
        id=AssetId("neale_lab_cfs_eur_raw"),
        trait="ME_CFS",
        project="neale_lab",
        sub_dir="raw",
        project_path=PurePath("20002_1482.gwas.imputed_v3.both_sexes.tsv.bgz"),
        read_spec=DataFrameReadSpec(format=DataFrameTextFormat(separator="\t")),
    ),
    url="https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/20002_1482.gwas.imputed_v3.both_sexes.tsv.bgz",
    md5_hash="b4d22a85b3a27265d39fb24ccb428461",
)
