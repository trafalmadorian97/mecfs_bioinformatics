from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import \
    KEATON_DBP_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_labeled_lead_variant import \
    SYNCOPE_LABELED_LEAD_VARIANT
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.analysis.syncope_s_ldsc import SYNCOPE_S_LDSC_TASKS
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_explode_and_scale import \
    SYNCOPE_SUMSTATS_EXPLODE_AND_SCALE
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_liftover import \
    AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS


def go():

    DEFAULT_RUNNER.run(
    SYNCOPE_S_LDSC_TASKS.get_terminal_tasks()
        ,

        must_rebuild_transitive=[
        ],
        incremental_save=True
    )


if __name__ == "__main__":
    go()
