from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_manhattan import MI_EUR_MANHATTAN
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import \
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.myocardial_infarction.raw.raw_mi_data import MILLION_VETERAN_MI_EUR_DATA_RAW
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import GENE_THESAURUS
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_cepo_specificity_matrix import \
    YU_DRG_CEPO_SPECIFICITY_MATRIX
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_counts_long import YU_DRG_COUNTS_LONG
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_counts_long_with_cell_type import \
    YU_DRG_COUNTS_LONG_WITH_CELL_TYPE
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_extracted_counts_data import \
    YU_DRG_EXTRACTED_COUNTS
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.processed.yu_drg_frac_specificity_matrix import \
    YU_DRG_FRAC_SPECIFICITY_MATRIX
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.raw.yu_drg_metadata_table import YU_DRG_METADATA_TABLE
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr15_54_925_638_locus_plalindindromes import \
#     DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr17_50_237_377_locus_palindromes import \
#     DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr1_174_128_548_locus_palindromes import \
#     DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr20_47_653_230_locus_palindromes import \
#     DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_26_215_000_locus_palindromes import \
#     DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP_PALINDROMES
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_97_505_620_locus_palindromes import \
#     DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES
#
# from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.without_palindromes.susie_finemap_decode_me_37_chr6_97_505_620_locus import \
#     DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP
from mecfs_bio.assets.reference_data.rna_seq_data.yu_drg.raw.yu_drg_raw_counts_rdata import YU_DRG_SRC1_RDATA


def run_mi_analysis():
    DEFAULT_RUNNER.run(
        # [MILLION_VETERAN_MI_EUR_DATA_RAW]+MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.get_terminal_tasks(),
        # [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task],
        [MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.heritability_task],
        # [MI_EUR_MANHATTAN],
        incremental_save=True,
        # must_rebuild_transitive=[MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS.manhattan_task]


    )


if __name__ == "__main__":
    run_mi_analysis()
