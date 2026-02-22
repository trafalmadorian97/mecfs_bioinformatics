"""
Script to run fine mapping analysis on DecodeME data.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr1_174_128_548_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_26_215_000_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_97_505_620_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr15_54_925_638_locus_plalindindromes import (
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr17_50_237_377_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr20_47_653_230_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES,
)


def run_fine_mapping_decode_me_analysis():
    """
    Function to fine map the main DecodeME loci using SUSIE.
    As an LD reference, use a UK biobank LD reference matrix from the Broad institute.
    """
    DEFAULT_RUNNER.run(
        (
            DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.terminal_tasks()
            + DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP_PALINDROMES.terminal_tasks()
            + DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.terminal_tasks()
            + DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES.terminal_tasks()
            + DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.terminal_tasks()
            + DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.terminal_tasks()
        ),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_fine_mapping_decode_me_analysis()
