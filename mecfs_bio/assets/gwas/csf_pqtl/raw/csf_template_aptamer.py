"""
The template aptamer for the CSF pQTL variant index: one arbitrary aptamer's
GWAS-SSF summary statistics, downloaded and converted to parquet.

ConstructCsfVariantIndexTask is templated off a single aptamer's variant set (the
CSF PLINK2 run is joint across aptamers, so one aptamer is representative). This
plays the role STACK_UKBBPPP_RABGAP1L plays for the PPP index.

The download is the raw .tsv.gz; PipeDataFrameTask reads the whole 7.3M-row file
with the polars backend (required for text sources) and rewrites it as parquet so
the index task can scan it cheaply. The read_spec on the download meta declares the
tab separator PipeDataFrameTask reads through.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.csf_database.gwas_catalog_url import (
    gwas_catalog_sumstats_url,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.constants.csf_database_constants import GcstAccession

# One arbitrary CSF aptamer (X13681.173 / CSNK2A2), with its published md5.
_TEMPLATE_ACCESSION = GcstAccession("GCST90421540")
_TEMPLATE_MD5 = "047befd46b553da2bcecf7c8faa91749"

CSF_TEMPLATE_APTAMER_DOWNLOAD = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="csf_pqtl",
        sub_group="western_2024",
        sub_folder=PurePath("raw"),
        id=AssetId("csf_template_aptamer_sumstats"),
        filename="csf_template_aptamer",
        extension=".tsv.gz",
        read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
    ),
    url=gwas_catalog_sumstats_url(_TEMPLATE_ACCESSION),
    md5_hash=_TEMPLATE_MD5,
)

CSF_TEMPLATE_APTAMER = PipeDataFrameTask.create(
    source_task=CSF_TEMPLATE_APTAMER_DOWNLOAD,
    asset_id="csf_template_aptamer_parquet",
    out_format=ParquetOutFormat(),
    pipes=[],
    backend="polars",
)
