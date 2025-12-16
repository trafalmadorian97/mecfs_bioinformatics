from mecfs_bio.assets.reference_data.db_snp.db_snp_build_37_as_parquet_unnest_ref import (
    PARQUET_DBSNP_37_UNNESTED,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe

PARQUET_DBSNP_37_UNNESTED_SORTED = PipeDataFrameTask.create(
    source_task=PARQUET_DBSNP_37_UNNESTED,
    asset_id="db_snp_build_37_parquet_unnested_sorted",
    out_format=ParquetOutFormat(),
    pipes=[SortPipe(by=["int_chrom", "POS"])],
    backend="ibis",
)
