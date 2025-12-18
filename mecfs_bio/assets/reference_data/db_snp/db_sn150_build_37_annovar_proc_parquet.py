from mecfs_bio.assets.reference_data.db_snp.db_snp150_build_37_annovar_proc import (
    DB_SNP150_ANNOVAR_PROC,
)
from mecfs_bio.build_system.task.compressed_csv_to_parquet_task import (
    CompressedCSVToParquetTask,
)

PARQUET_DBSNP150_37_ANNOVAR_PROC = CompressedCSVToParquetTask.create(
    csv_task=DB_SNP150_ANNOVAR_PROC,
    asset_id="db_snp150_annovar_proc_parquet",
    source_compression=None,
    type_dict={"CHROM": "VARCHAR"},
)
