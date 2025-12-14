from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.red_blood_cell_volume.million_veterans.analysis.rbc_volume_magma_task_generator import \
    MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS


def run_red_blood_cell_volume_analysis():
    """
    Analyze GWAS of red blood cell volume from million veterans
    Includes:
    - Apply MAGMA to estimate significant genes
    - Combine Magma results with GTEx to estimate signifiant tissues


    """
    DEFAULT_RUNNER.run(
        [MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS.inner.bar_plot_task,
         MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS.inner.gene_analysis_task],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_red_blood_cell_volume_analysis()
