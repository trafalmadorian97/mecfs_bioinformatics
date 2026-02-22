"""
Task to Assign rsids to SNPs using the ucsc hg19 sbSNP build 151 reference file.

Note that this approach only allows assignment of rsids to SNPs.  Non-SNP variants are dropped.

"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_37_parquet_file import (
    DECODE_ME_GWAS_1_LIFTOVER_TO_37_PARQUET,
)
from mecfs_bio.assets.reference_data.db_snp.snp151_build37_parquet import (
    GENOME_ANNOTATION_DATABASE_BUILD_37_PARQUET,
)
from mecfs_bio.build_system.reference.schemas.chrom_rename_rules import (
    CHROM_RENAME_RULES,
)
from mecfs_bio.build_system.reference.schemas.hg19_snp151_schema_valid_choms import (
    HG19_SNP151_VALID_CHROMS,
)
from mecfs_bio.build_system.task.assign_rsids_via_snp151_task import (
    AssignRSIDSToSNPsViaSNP151Task,
)

DECODE_ME_GWAS_1_LIFTOVER_TO_37_WITH_RSID = AssignRSIDSToSNPsViaSNP151Task.create(
    snp151_database_file_task=GENOME_ANNOTATION_DATABASE_BUILD_37_PARQUET,
    raw_snp_data_task=DECODE_ME_GWAS_1_LIFTOVER_TO_37_PARQUET,
    asset_id="decode_me_gwas_1_liftover_to_build_37_with_rsid",
    valid_chroms=HG19_SNP151_VALID_CHROMS,
    chrom_replace_rules=CHROM_RENAME_RULES,
)
