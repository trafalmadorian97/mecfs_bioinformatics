"""
Rough experimental scripts pertaining to analysis of ME/CFS data
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.multi_gwas.genetic_correlation.ct_ldsc_plot import CT_LDSC_INITIAL_PLOT


def run_initial_decode_me_analysis():
    DEFAULT_RUNNER.run(
        # [
    # [DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP.susie_finemap_2_credible_set_plot]
    #     DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP.terminal_tasks()
    # DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP.terminal_tasks()
   # DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP.terminal_tasks()

        #     DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP.terminal_tasks(),
        (
            # [YU_DRG_SRC1_RDATA]
            [
                CT_LDSC_INITIAL_PLOT
                # YU_DRG_FRAC_SPECIFICITY_MATRIX
                # CT_LDSC_INITIAL_PLOT
                # CT_LDSC_INITIAL
                # JOHNSON_DRG_MAGMA_CEPO_BAR_PLOT,
                # JOHNSON_DRG_MAGMA_FRAC_BAR_PLOT
                # JOHNSON_DRG_MAGMA_FRAC_BAR_PLOT
                # YU_DRG_CEPO_SPECIFICITY_MATRIX
                # YU_DRG_EXTRACTED_COUNTS,
             # YU_DRG_METADATA_TABLE,
             # GENE_THESAURUS,
             # YU_DRG_COUNTS_LONG,
             # YU_DRG_COUNTS_LONG_WITH_CELL_TYPE
             ]
        # DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES.terminal_tasks()
        #     DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP_PALINDROMES.terminal_tasks()
    # DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.terminal_tasks()
    #     DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES.terminal_tasks()+
    #     DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES.terminal_tasks()
    #     DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES.terminal_tasks()
    #     DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES.terminal_tasks()
    ),

    # DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROME.terminal_tasks(),
        # MILLION_VETERAN_LDL_EUR_DATA_RAW
            # DECODE_ME_RABGAP1L_REGION_PLOT_37,
            # DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task
            # CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES
            # DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS,
            # DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS_STACKPLOT
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES,
            # CHR1_173000001_17600000_UKBB_LD_MATRIX_DOWNLOAD,
            # CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD
         # DECODE_ME_BTN1A1_REGION_PLOT_37
         #    DECODE_ME_BTN1A1_REGION_PLOT_37
         #    ROADMAP_CELL_TYPE_CATEGORIES_FOR_LDSC

         # ],
        # DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks(),
        # DECODE_ME_BASIC_CIS_PQTL_MR.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            CT_LDSC_INITIAL_PLOT.genetic_corr_source.task
            # CT_LDSC_INITIAL,
            # CT_LDSC_INITIAL_PLOT
            # CT_LDSC_INITIAL
            # YU_DRG_CEPO_SPECIFICITY_MATRIX
            # YU_DRG_FRAC_SPECIFICITY_MATRIX

            # YU_DRG_CEPO_SPECIFICITY_MATRIX
        ]
        # must_rebuild_transitive= DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP.terminal_tasks(),
            # DECODE_ME_GWAS_37_CHR6_26_215_000_FINEMAP.harmonized_sumstats_task
            # DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS_STACKPLOT
            # DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD
            # DECODE_ME_BASIC_CIS_PQTL_MR.multiple_testing_task
            # DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
            # DECODE_ME_MASTER_GENE_LIST_WITH_GGET
            # DECODE_ME_HBA_MAGMA_TASKS.magma_independent_cluster_plot

    )


if __name__ == "__main__":
    run_initial_decode_me_analysis()
