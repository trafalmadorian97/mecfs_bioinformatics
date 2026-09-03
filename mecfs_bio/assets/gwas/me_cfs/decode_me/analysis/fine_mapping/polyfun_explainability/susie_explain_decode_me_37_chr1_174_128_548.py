"""Demonstrator: polyfun-vs-uniform explainability at the DecodeME chr1:174.1Mb
locus.

Wires the polyfun explainability outer generator at the same locus, sumstats,
sample size, palindrome strategy, and chrom range as the existing
with_palindromes_and_polyfun_precomputed_prior chr1_174 fine-mapping module, so
the eight SUSIE runs here operate on the same harmonized inputs. The result is
an 8-run outer group (four run configs, each a matched uniform/polyfun pair)
plus a contrast and plot task per pair.
"""

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    generate_assets_polyfun_explain_fine_map,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    ChromRange,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

POLYFUN_EXPLAIN_CHR1_174 = generate_assets_polyfun_explain_fine_map(
    chrom=1,
    pos=174_128_548,
    build_37_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    base_name="decode_me_polyfun_explain",
    sumstats_pipe=IdentityPipe(),
    sample_size_or_effect_sample_size=int(
        4 / (1 / 15_579 + 1 / 259_909)
    ),  # 4/(1/cases + 1/controls)
    palindrome_strategy="keep",
    chrom_range=ChromRange(1, 173_500_000, 174_500_000),
)
