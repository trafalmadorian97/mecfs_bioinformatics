from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr15_54_925_638_locus_plalindindromes import \
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr17_50_237_377_locus_palindromes import \
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr1_174_128_548_locus_palindromes import \
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr20_47_653_230_locus_palindromes import \
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_97_505_620_locus_palindromes import \
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES
from mecfs_bio.figures.key_scripts.regenerate_figures import regenerate_figures


def go():
    regenerate_figures(
        [

            DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.susie_base_credible_set_parquet_table,
            DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.susie_base_credible_set_parquet_table,
            DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.susie_2_credible_set_parquet_table,
            DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES.susie_base_credible_set_parquet_table,
            DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.susie_base_credible_set_parquet_table,
            DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_base_credible_set_parquet_table,
            DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.susie_strict_credible_set_parquet_table,

        ]
    )

if __name__ == '__main__':
    go()