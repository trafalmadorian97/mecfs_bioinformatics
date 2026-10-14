from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_build_37_with_rsid import (
    DECODE_ME_GWAS_1_LIFTOVER_TO_37_WITH_RSID,
)
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask

DECODE_ME_GWAS_1_BUILD_37_MAGMA_SNP_LOCS = (
    MagmaSNPFileTask.create_for_magma_snp_pos_file(
        gwas_parquet_with_rsids_task=DECODE_ME_GWAS_1_LIFTOVER_TO_37_WITH_RSID,
        asset_id="decode_me_gwas_1_build_37_magma_snp_locs",
    )
)
