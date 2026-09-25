"""Polyfun-vs-uniform explainability at the DecodeME chr6:97.5Mb locus, with the
DecodeME-specific L2-regularized S-LDSC prior (PolyFun Approach 2) in place of the
precomputed PolyFun prior.

Same locus, sumstats, sample size, palindrome strategy, and chrom range as
polyfun_explainability.susie_explain_decode_me_37_chr6_97_505_620, so the two sets of
credible sets differ only in the prior. The contrast attributes the prior to
annotations with the DecodeME tau that scored chr6 (fit on the odd
chromosomes), which explains the unfloored prior exactly.
"""

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    generate_assets_polyfun_explain_fine_map,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.polyfun_approach_2_l2_sldsc_prior.decode_me_l2_sldsc_prior import (
    DECODE_ME_L2_SLDSC_PRIOR_SOURCE,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    ChromRange,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)

POLYFUN_EXPLAIN_CHR6_97_L2_SLDSC_PRIOR = generate_assets_polyfun_explain_fine_map(
    chrom=6,
    pos=97_505_620,
    build_37_sumstats_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    base_name="decode_me_polyfun_explain_l2_sldsc_prior",
    sumstats_pipe=IdentityPipe(),
    sample_size_or_effect_sample_size=DECODE_ME_PREVALENCE_INFO.effective_sample_size,
    palindrome_strategy="keep",
    chrom_range=ChromRange(6, 97_500_000, 99_000_000),
    secondary_position_from_snpid=SecondaryPositionFromSnpid(build_label="hg38"),
    prior_source=DECODE_ME_L2_SLDSC_PRIOR_SOURCE,
)
