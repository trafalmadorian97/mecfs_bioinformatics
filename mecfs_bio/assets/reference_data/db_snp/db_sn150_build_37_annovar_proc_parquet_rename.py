"""
Task to rename the chromosomes in the annovar version of dbsnp 150 to match
those used by gwaslab
"""

from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RD,
)
from mecfs_bio.build_system.reference.schemas.chrom_rename_rules import (
    CHROM_RENAME_DF_NUMERIC_X_Y_M,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.join_with_memory_table_pipe import (
    JoinWithMemTablePipe,
)
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe

PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME = PipeDataFrameTask.create(
    asset_id="db_snp150_annovar_proc_parquet_rename",
    source_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RD,
    pipes=[
        JoinWithMemTablePipe(
            CHROM_RENAME_DF_NUMERIC_X_Y_M,
            keys_left=["CHROM"],
            keys_right=["CHROM"],
            backend="ibis",
        ),
        SortPipe(by=["int_chrom", "POS"]),
    ],
    out_format=ParquetOutFormat(),
    backend="ibis",
)

PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_RD = DiscardDepsWrapper(
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME
)
