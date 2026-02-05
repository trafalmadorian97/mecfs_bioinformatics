
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_cis_pqtl_mr import DECODE_ME_BASIC_CIS_PQTL_MR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_region_plot_BTN1A1_locus_37 import \
    DECODE_ME_BTN1A1_REGION_PLOT_37
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_region_plot_rabgap1l_locus_37 import \
    DECODE_ME_RABGAP1L_REGION_PLOT_37
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.susie_finemap_decode_me_37_chr1_173_locus import \
    DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import \
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.prep_for_fine_mapping.chr1_173_locus.harmonize_with_polyfun_reference_alleles import \
    DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.prep_for_fine_mapping.chr1_173_locus.harmonize_with_polyfun_reference_rsid import \
    DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.custom.roadmap_cell_type_categorization import \
    ROADMAP_CELL_TYPE_CATEGORIES_FOR_LDSC
from mecfs_bio.assets.reference_data.ukbb_ld_matrices.from_polyfun.chr1_173000001.chr1_173000001_176000001_labels import \
    CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD


def run_initial_decode_me_analysis():
    DEFAULT_RUNNER.run(
        [
            # DECODE_ME_RABGAP1L_REGION_PLOT_37,
            # DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task
            # CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES
            DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS
         # DECODE_ME_BTN1A1_REGION_PLOT_37
         #    DECODE_ME_BTN1A1_REGION_PLOT_37
         #    ROADMAP_CELL_TYPE_CATEGORIES_FOR_LDSC

         ],
        # DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks(),
        # DECODE_ME_BASIC_CIS_PQTL_MR.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[

            DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS
            # DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD
            # DECODE_ME_BASIC_CIS_PQTL_MR.multiple_testing_task
            # DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
            # DECODE_ME_MASTER_GENE_LIST_WITH_GGET
            # DECODE_ME_HBA_MAGMA_TASKS.magma_independent_cluster_plot
        ],
    )


if __name__ == "__main__":
    run_initial_decode_me_analysis()
