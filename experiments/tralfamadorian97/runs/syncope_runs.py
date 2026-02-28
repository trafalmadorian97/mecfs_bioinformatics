from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import \
    KEATON_DBP_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_hba_magma_task_generator import \
    HBA_MAGMA_AEGISDOTTIR_SYCNOPE_GWAS
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_labeled_lead_variant import \
    SYNCOPE_LABELED_LEAD_VARIANT
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_s_ldsc import SYNCOPE_S_LDSC_TASKS
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_explode_and_scale import \
    SYNCOPE_SUMSTATS_EXPLODE_AND_SCALE
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_liftover import \
    AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS


from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.sumstats_liftover_dump_to_parquet import AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET
def go():

    DEFAULT_RUNNER.run(
        # [AEGISDOTTIR_SUMSTATS_DUMP_TO_PARQUET]
    HBA_MAGMA_AEGISDOTTIR_SYCNOPE_GWAS.terminal_tasks()
    # SYNCOPE_S_LDSC_TASKS.get_terminal_tasks()
        ,

        must_rebuild_transitive=[
        ],
        incremental_save=True
    )


if __name__ == "__main__":
    go()
