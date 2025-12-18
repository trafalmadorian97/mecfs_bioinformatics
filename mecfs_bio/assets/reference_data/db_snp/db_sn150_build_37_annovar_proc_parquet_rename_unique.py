from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe

PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE = PipeDataFrameTask.create(
    asset_id="db_snp150_annovar_proc_parquet_rename_unique",
    source_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME,
    pipes=[UniquePipe(by=["int_chrom", "POS", "ALT", "REF"])],
    out_format=ParquetOutFormat(),
    backend="ibis",
)
