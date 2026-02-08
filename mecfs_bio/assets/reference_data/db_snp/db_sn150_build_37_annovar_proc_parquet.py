from mecfs_bio.assets.reference_data.db_snp.db_snp150_build_37_annovar_proc import (
    DB_SNP150_ANNOVAR_PROC,
    DB_SNP150_ANNOVAR_PROC_RD,
)
from mecfs_bio.build_system.task.compressed_csv_to_parquet_task import (
    CompressedCSVToParquetTask,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper

PARQUET_DBSNP150_37_ANNOVAR_PROC = CompressedCSVToParquetTask.create(
    csv_task=DB_SNP150_ANNOVAR_PROC_RD,
    asset_id="db_snp150_annovar_proc_parquet",
    source_compression=None,
    type_dict={"CHROM": "VARCHAR"},
)


PARQUET_DBSNP150_37_ANNOVAR_PROC_NON_RD = CompressedCSVToParquetTask.create(
    csv_task=DB_SNP150_ANNOVAR_PROC,
    asset_id="db_snp150_annovar_proc_parquet",
    source_compression=None,
    type_dict={"CHROM": "VARCHAR"},
)

PARQUET_DBSNP150_37_ANNOVAR_PROC_RD = DiscardDepsWrapper(
    PARQUET_DBSNP150_37_ANNOVAR_PROC
)
