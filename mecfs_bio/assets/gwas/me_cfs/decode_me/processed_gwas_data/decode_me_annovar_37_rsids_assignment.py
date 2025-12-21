from mecfs_bio.asset_generator.annovar_37_basic_rsid_assignment import (
    annovar_37_basic_rsid_assignment,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import (
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
)

DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED = annovar_37_basic_rsid_assignment(
    sumstats_task=DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37, base_name="decode_me_gwas_1"
)
