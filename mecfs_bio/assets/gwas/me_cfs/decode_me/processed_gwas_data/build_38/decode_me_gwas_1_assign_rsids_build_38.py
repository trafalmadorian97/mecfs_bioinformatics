from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_minimal_processing import (
    DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
)
from mecfs_bio.assets.reference_data.db_snp.db_snp_vcf_file_with_index_build_38 import (
    DB_SNP_VCF_FILE_WITH_INDEX_BUILD_38_DIR,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_assign_rsids_via_sweep_task import (
    GWASLabAssignRSIDSViaSweepTask,
)

DECODE_ME_GWAS_1_ASSIGN_RSIDS_BUILD_38 = GWASLabAssignRSIDSViaSweepTask.create(
    id="decode_me_gwas_1_assign_rsids_build_38",
    sumstats_task=DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
    vcf_dir_task=DB_SNP_VCF_FILE_WITH_INDEX_BUILD_38_DIR,
    vcf_filename="GCF_000001405.40.gz",
)
