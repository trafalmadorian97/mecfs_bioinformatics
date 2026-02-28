from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import \
    KEATON_DBP_STANDARD_ANALYSIS


def go():

    DEFAULT_RUNNER.run(
        [
            KEATON_DBP_STANDARD_ANALYSIS.labeled_lead_variant_tasks.raw_sumstats_task

        ],

        must_rebuild_transitive=[
        ],
    )


if __name__ == "__main__":
    go()
