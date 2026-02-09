"""
Task to deduplicate variants in the annovar version of dbSNP150
This is necessary because dbSNP can contain multiple rsids that refer to the same genetic variant.

The current deduplication strategy is just to take the largest rsid in the string ordering sense
This is not optimal.

The optimal choice is to select the rsid that will match what is expected by the downstream task.
For example, if we are going to perform S-LDSC, for each variant we should choose the rsid that matches the one in
the S-LDSC reference data.

"""

from pathlib import PurePath

from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_NON_RD,
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_RD,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.duckdb_mem_limit_pipe import DuckdbMemLimitPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe

PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE = DiscardDepsWrapper(
    PipeDataFrameTask.create(
        asset_id="db_snp150_annovar_proc_parquet_rename_unique",
        source_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_RD,
        pipes=[
            DuckdbMemLimitPipe(limit_gb=4),
            UniquePipe(
                by=["int_chrom", "POS", "ALT", "REF"],
                keep="last",
                order_by=["int_chrom", "POS", "ALT", "REF", "rsid"],
            ),
        ],
        out_format=ParquetOutFormat(),
        backend="duckdb",
    )
)


PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_DIRECT_DOWNLOAD = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="db_snp_reference_data",
        sub_group="build_37",
        sub_folder=PurePath("annovar"),
        id="db_snp150_annovar_proc_parquet_unique_direct_download",
        filename="db_snp150_annovar_proc_parquet_unique_direct_download",
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    ),
    url="https://www.dropbox.com/scl/fi/6oglew7b7eh3zk0bn9gen/db_snp150_annovar_proc_parquet_rename_unique.parquet?rlkey=geuwlsha5icnthu3bxoyvchnm&dl=1",
    md5_hash="1e990169a3ff48bae01f53cc69336be3",
)


PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_NON_RD = PipeDataFrameTask.create(
    asset_id="db_snp150_annovar_proc_parquet_rename_unique",
    source_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_NON_RD,
    pipes=[
        DuckdbMemLimitPipe(limit_gb=4),
        UniquePipe(
            by=["int_chrom", "POS", "ALT", "REF"],
            keep="last",
            order_by=["int_chrom", "POS", "ALT", "REF", "rsid"],
        ),
    ],
    out_format=ParquetOutFormat(),
    backend="duckdb",
)
