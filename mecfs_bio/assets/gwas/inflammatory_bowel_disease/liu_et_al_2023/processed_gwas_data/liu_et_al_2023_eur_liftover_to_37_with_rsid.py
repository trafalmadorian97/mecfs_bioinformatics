from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_parquet import (
    LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_PARQUET,
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

LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_WITH_RSID = AssignRSIDSToSNPsViaSNP151Task.create(
    snp151_database_file_task=GENOME_ANNOTATION_DATABASE_BUILD_37_PARQUET,
    raw_snp_data_task=LIU_ET_AL_2023_IBD_LIFTOVER_TO_37_PARQUET,
    asset_id="liu_et_al_2023_ibd_eur_liftover_to_37_parquet_with_rsid",
    valid_chroms=HG19_SNP151_VALID_CHROMS,
    chrom_replace_rules=CHROM_RENAME_RULES,
)
