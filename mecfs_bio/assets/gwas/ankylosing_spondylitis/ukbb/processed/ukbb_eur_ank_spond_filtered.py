"""
UKBB sumstats are from whole genome sequencing, and so contain a huge number of variants
Most downstream tools expect genotype data, so we can safely throw out many variants

Note that the data used as a filter is from million veterans and is build 37.  However, since rsids do not change when build changes, this does not matter
"""

from mecfs_bio.assets.gwas.ankylosing_spondylitis.million_veterans.processed.mv_eur_ank_spond_sumstats_dump_to_parquet import (
    MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET,
)
from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_eur_ank_spond_parquet import (
    UKBB_ANK_SPOND_PARQUET,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_RSID_COL,
)

FILTERED_UKBB_ANK_SPOND = JoinDataFramesTask.create_from_result_df(
    backend="ibis",
    asset_id="ukbb_ank_spond_sumstats_filtered",
    result_df_task=UKBB_ANK_SPOND_PARQUET,
    reference_df_task=MILLION_VETERANS_ANK_SPOND_SUMSTATS_37_DUMP_TO_PARQUET,
    how="inner",
    left_on=["rsid"],
    right_on=[GWASLAB_RSID_COL],
    out_format=ParquetOutFormat(),
    df_2_pipe=CompositePipe(
        [
            SelectColPipe(
                [GWASLAB_RSID_COL],
            ),
            UniquePipe([GWASLAB_RSID_COL], keep="any", order_by=[GWASLAB_RSID_COL]),
        ]
    ),
)
