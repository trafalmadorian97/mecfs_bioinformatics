from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_multisite_pain import \
    JOHNSTON_ET_AL_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.magma.johnston_gene_level_magma_plot import \
    JOHNSTON_MAGMA_GENE_PLOT
from mecfs_bio.build_system.scheduler.topological_scheduler import TopologicalSchedulerSettings


def run_pain():
    DEFAULT_RUNNER.run(


            # JOHNSTON_MAGMA_GENE_PLOT
            JOHNSTON_ET_AL_PPP_RG_CIS_EXCLUDED.get_terminal_tasks()
        ,


        incremental_save=True,
        must_rebuild_transitive=[
        ],
        settings=TopologicalSchedulerSettings(
            print_progress=False
        )


    )


if __name__ == "__main__":
    run_pain()
