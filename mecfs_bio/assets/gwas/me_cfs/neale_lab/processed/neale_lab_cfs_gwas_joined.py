from mecfs_bio.assets.gwas.me_cfs.neale_lab.raw.neale_lab_cfs_gwas_download import (
    NEALE_LAB_CFS_RAW,
)
from mecfs_bio.assets.reference_data.neale_lab_reference_data.neale_lab_variant_list import (
    NEALE_LAB_VARIANTS_REFERENCE,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

NEALE_LAB_CFS_JOINED = JoinDataFramesTask.create_from_result_df(
    asset_id="neale_lab_cfs_gwas_joined",
    result_df_task=NEALE_LAB_CFS_RAW,
    reference_df_task=NEALE_LAB_VARIANTS_REFERENCE,
    how="inner",
    left_on=["variant"],
    right_on=["variant"],
    out_format=ParquetOutFormat(),
    df_1_pipe=CompositePipe(
        [FilterRowsByValue("low_confidence_variant", valid_values=[False])],
    ),
    df_2_pipe=CompositePipe(
        [
            SelectColPipe(
                [
                    "variant",
                    "chr",
                    "pos",
                    "ref",
                    "alt",
                    "rsid",
                    "consequence",
                    "consequence_category",
                    "info",
                ]
            )
        ]
    ),
)
