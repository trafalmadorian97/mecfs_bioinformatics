from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.varicella_zoster_antibodies.butler_laporte_et_al.analysis.butler_laporte_vzv_standard_analysis import \
    BUTLER_LAPORTE_ET_AL_VZV_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.varicella_zoster_antibodies.butler_laporte_et_al.raw.vzv_igg_raw import \
    BUTLER_LAPORTE_ET_AL_VZV_SUMSTATS


def go():

    DEFAULT_RUNNER.run(
        BUTLER_LAPORTE_ET_AL_VZV_STANDARD_ANALYSIS.get_terminal_tasks()
        # [BUTLER_LAPORTE_ET_AL_VZV_SUMSTATS]
        ,
        must_rebuild_transitive=[
        ],
        incremental_save=True,

    )


if __name__ == "__main__":
    go()
