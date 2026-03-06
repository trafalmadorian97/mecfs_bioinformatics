from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.raw.raw_ukbb_eur_ank_spond_data import (
    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

UKBB_ANK_SPOND_PARQUET = PipeDataFrameTask.create(
    source_task=UK_BIOBANK_ANKYLOSING_SPONDYLITIS_EUR_DATA_RAW,
    asset_id="ubb_ank_spond_parquet",
    out_format=ParquetOutFormat(),
    pipes=[IdentityPipe()],
    backend="polars",
)
