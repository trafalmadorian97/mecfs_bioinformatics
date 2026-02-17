from mecfs_bio.asset_generator.fine_mapping_asset_generator import (
    generate_assets_broad_ukbb_fine_map,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    ChromRange,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES = generate_assets_broad_ukbb_fine_map(
    chrom=17,
    pos=50_237_377,
    build_37_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    base_name="decode_me",
    sumstats_pipe=IdentityPipe(),
    sample_size_or_effect_sample_size=int(
        4 / (1 / 15_579 + 1 / 259_909)
    ),  # 4/(1/cases + 1/controls)
    palindrome_strategy="keep",
    chrom_range=ChromRange(17, 50_000_000, 51_000_000),
)
