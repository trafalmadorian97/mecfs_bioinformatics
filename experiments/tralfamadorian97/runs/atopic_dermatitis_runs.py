
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.atopic_dermatitis.budu_aggrey_et_al_2023.analysis.budu_aggrey_analysis import \
    BUDU_AGGREY_ET_AL_ATOPIC_DERMATITIS_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.atopic_dermatitis.budu_aggrey_et_al_2023.raw.raw_budu_aggrey_atopic_dermatitis_data import \
    BUDU_AGGREY_ET_AL_ATOPIC_DERMATITIS_RAW


def run_initial_atopic_derm_analysis():
    DEFAULT_RUNNER.run(
        (
            BUDU_AGGREY_ET_AL_ATOPIC_DERMATITIS_STANDARD_ANALYSIS.get_terminal_tasks()

    ),

        incremental_save=True,
        must_rebuild_transitive=[
        ]

    )


if __name__ == "__main__":
    run_initial_atopic_derm_analysis()
