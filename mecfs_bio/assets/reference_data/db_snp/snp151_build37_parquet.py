from mecfs_bio.assets.reference_data.db_snp.snp151_build_37 import (
    GENOME_ANNOTATION_DATABASE_BUILD_37,
)
from mecfs_bio.build_system.task.compressed_csv_to_parquet_task import (
    CompressedCSVToParquetTask,
)

GENOME_ANNOTATION_DATABASE_BUILD_37_PARQUET = CompressedCSVToParquetTask.create(
    csv_task=GENOME_ANNOTATION_DATABASE_BUILD_37,
    asset_id="genome_annotation_database_build_37_as_parquet",
    source_compression="gzip",
)
