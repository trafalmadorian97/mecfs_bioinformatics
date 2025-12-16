from mecfs_bio.assets.reference_data.db_snp.db_snp_build_37_as_parquet import (
    PARQUET_DBSNP_37,
)
from mecfs_bio.build_system.reference.nc_chrom_mapping import (
    NC_CHROM_MAPPING,
    NC_CHROM_TABLE,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.pipes.join_with_memory_table_pipe import (
    JointWithMemTablePipe,
)
from mecfs_bio.build_system.task.pipes.str_split_col import SplitColPipe
from mecfs_bio.build_system.task.pipes.unnest_pipe import UNNestPipe

PARQUET_DBSNP_37_UNNESTED = PipeDataFrameTask.create(
    source_task=PARQUET_DBSNP_37,
    asset_id="db_snp_build_37_parquet_unnested",
    out_format=ParquetOutFormat(),
    pipes=[
        SplitColPipe(col_to_split="ALT", split_by=",", new_col_name="ALT"),
        UNNestPipe(col_to_unnest="ALT"),
        FilterRowsByValue(
            target_column="CHROM",
            valid_values=list(NC_CHROM_MAPPING.keys()),
        ),
        JointWithMemTablePipe(
            mem_table=NC_CHROM_TABLE,
            keys_left=["CHROM"],
            keys_right=["nc_chrom"],
        ),
    ],
    backend="polars",
)
