"""
Asset generator for assigning RSIDS to genome build-37 GWS datasets.
"""

import narwhals
from attrs import frozen

from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename_unique import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GwasLabTransformSpec,
    GWASLabVCFRef,
    HarmonizationOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import (
    GWASLabTransformSumstatsTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe


@frozen
class RSIDAssignmentTaskGroup:
    harmonize: Task
    dump_parquet: Task
    join_task: Task

def annovar_37_basic_rsid_assignment(
    sumstats_task: Task,
    base_name: str,
    use_gwaslab_rsids_convention:bool=False
) -> RSIDAssignmentTaskGroup:
    """
    Asset generator that creates a chain of tasks to assign rsids to existing build 37 sumstats datasets using the annovar dbSNP reference data
    """
    harmonized_task = GWASLabTransformSumstatsTask.create_from_source_task(
        sumstats_task,
        asset_id=base_name + "__harmonized",
        spec=GwasLabTransformSpec(
            harmonize_options=HarmonizationOptions(
                ref_infer=GWASLabVCFRef(name="1kg_eur_hg19", ref_alt_freq="AF"),
                ref_seq="ucsc_genome_hg19",
                check_ref_files=True,
                drop_missing_from_ref=True,
                cores=4,
            )
        ),
    )
    dump_parquet_task = GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=harmonized_task,
        asset_id=base_name + "_harmonized_dump_to_parquet",
        sub_dir="processed",
    )
    if use_gwaslab_rsids_convention:
        out_pipe = RenameColPipe(
            old_name="rsid",
            new_name="rsID"
        )
    else:
        out_pipe = IdentityPipe()
    join_with_rsid_task = JoinDataFramesTask.create_from_result_df(
        asset_id=base_name + "_assign_rsids_via_dbsnp150",
        result_df_task=dump_parquet_task,
        reference_df_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE,
        left_on=["CHR", "POS", "EA", "NEA"],
        right_on=["int_chrom", "POS", "ALT", "REF"],
        out_format=ParquetOutFormat(),
        how="inner",
        df_1_pipe=CompositePipe(
            [
                CastPipe(
                    target_column="EA",
                    type=narwhals.dtypes.String(),
                    new_col_name="EA",
                ),
                CastPipe(
                    target_column="NEA",
                    type=narwhals.dtypes.String(),
                    new_col_name="NEA",
                ),
            ]
        ),
        backend="ibis",
        out_pipe=out_pipe

    )
    group =RSIDAssignmentTaskGroup(
        harmonize=harmonized_task,
        dump_parquet=dump_parquet_task,
        join_task=join_with_rsid_task,
    )
    return group
