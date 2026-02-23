
import narwhals.dtypes

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar import \
    LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_dump_to_parquet import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
)
from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename_unique import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_NON_RD,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabCreateSumstatsTask
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe


LIU_ET_AL_2023_RSIDS_FROM_ANNOVAR_SUMSATS=GWASLabCreateSumstatsTask(
    df_source_task=LIU_ET_AL_2023_ASSIGN_RSID_VIA_SNP150_ANNOVAR,
    asset_id=AssetId("liu_et_al_2023_rsids_from_annovar_sumstats"),
    basic_check=True,
    genome_build="19",
    fmt="gwaslab"
)