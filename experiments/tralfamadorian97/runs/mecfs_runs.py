
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_cis_pqtl_mr import DECODE_ME_BASIC_CIS_PQTL_MR
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_region_plot_BTN1A1_locus_37 import \
    DECODE_ME_BTN1A1_REGION_PLOT_37
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_region_plot_rabgap1l_locus_37 import \
    DECODE_ME_RABGAP1L_REGION_PLOT_37


def run_initial_decode_me_analysis():
    DEFAULT_RUNNER.run(
        [
            # DECODE_ME_RABGAP1L_REGION_PLOT_37,

         # DECODE_ME_BTN1A1_REGION_PLOT_37
            DECODE_ME_BTN1A1_REGION_PLOT_37
         ],
        # DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks(),
        # DECODE_ME_BASIC_CIS_PQTL_MR.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            # DECODE_ME_BASIC_CIS_PQTL_MR.multiple_testing_task
            # DECODE_ME_MASTER_GENE_LIST_AS_MARKDOWN,
            # DECODE_ME_MASTER_GENE_LIST_WITH_GGET
            # DECODE_ME_HBA_MAGMA_TASKS.magma_independent_cluster_plot
        ],
    )


if __name__ == "__main__":
    run_initial_decode_me_analysis()
