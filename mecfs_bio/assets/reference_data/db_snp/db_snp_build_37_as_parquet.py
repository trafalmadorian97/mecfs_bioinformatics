from mecfs_bio.assets.reference_data.db_snp.db_snp_vcf_file_for_rsid_assignment_build_37 import (
    DB_SNP_VCF_FILE_BUILD_37,
)
from mecfs_bio.build_system.task.compressed_csv_to_parquet_task import (
    CompressedCSVToParquetTask,
)

PARQUET_DBSNP_37 = CompressedCSVToParquetTask.create(
    csv_task=DB_SNP_VCF_FILE_BUILD_37,
    asset_id="db_snp_vcf_build_37_parquet",
    source_compression="gzip",
    select_list=[
        "CHROM",
        "POS",
        "ID",
        "REF",
        "ALT",
    ],
)
