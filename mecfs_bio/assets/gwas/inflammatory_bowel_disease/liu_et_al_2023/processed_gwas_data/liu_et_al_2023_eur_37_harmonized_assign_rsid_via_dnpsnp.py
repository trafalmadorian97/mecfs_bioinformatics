import narwhals.dtypes

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_dump_to_parquet import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
)
from mecfs_bio.assets.reference_data.db_snp.db_snp_build_37_as_parquet_unnest_ref_sorted import (
    PARQUET_DBSNP_37_UNNESTED_SORTED,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe

LIU_ET_AL_2023_ASSIGN_RSID_VIA_DBSNP = JoinDataFramesTask.create_from_result_df(
    asset_id="liu_et_al_2023_ibd_eur_harmonize_assign_rsids_via_dnpsnp",
    result_df_task=LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
    reference_df_task=PARQUET_DBSNP_37_UNNESTED_SORTED,
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
)
