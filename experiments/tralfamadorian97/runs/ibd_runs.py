
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.cis_pqtl_mr import LIU_ET_AL_CIS_PQTL_MR
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.magma.liu_et_al_2023_eur_37_hba_magma_analysis import (
    IBD_HBA_MAGMA_TASKS,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.analysis.magma.liu_et_al_2023_eur_37_specific_tissue_bar_plot import (
    LIU_ET_AL_IBD_EUR_37_SPECIFIC_TISSUE_ANALYSIS_BAR_PLOT,
)


def run_ibd():
    DEFAULT_RUNNER.run(
        LIU_ET_AL_CIS_PQTL_MR.terminal_tasks(),

        # + LIU_ET_AL_S_LSDC_FROM_SNP_150.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
            LIU_ET_AL_CIS_PQTL_MR.md_task
        ]
        # must_rebuild_transitive=[
        #     LIU_ET_AL_CIS_PQTL_MR.mr_task
        # ]
        # must_rebuild_transitive=[IBD_HBA_MAGMA_TASKS.magma_independent_cluster_plot],
    )


if __name__ == "__main__":
    run_ibd()
