"""
Task to deduplicate variants in the annovar version of dbSNP150
This is necessary because dbSNP can contain multiple rsids that refer to the same genetic variant.

The current deduplication strategy is just to take the largest rsid in the string ordering sense
This is not optimal.

The optimal choice is to select the rsid that will match what is expected by the downstream task.
For example, if we are going to perform S-LDSC, for each variant we should choose the rsid that matches the one in
the S-LDSC reference data.

"""

from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_NON_RD,
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_RD,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
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
