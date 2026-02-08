"""
Task to use the version of dbSNP 150 downloaded from annovar to assign rsids to variants from the Liu et al. study.

Because the Annovar version of dbSNP is left normalized, it is more likely to contain matches
for the genetic variants in the Liu et al. study
"""

import narwhals.dtypes

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_dump_to_parquet import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
)
from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename_unique import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE, PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_NON_RD,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe

LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR = (
    JoinDataFramesTask.create_from_result_df(
        asset_id="liu_et_al_2023_ibd_eur_harmonize_assign_rsids_via_snp150_annovar",
        result_df_task=LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
        reference_df_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_NON_RD,
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
    )
)
