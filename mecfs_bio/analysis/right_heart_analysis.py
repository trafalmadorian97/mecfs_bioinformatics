from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.imaging_derived_heart_phenotypes.pirruccello_et_al_2022.processed.extracted_rvef_data import \
    PIRRUCCELLO_EXTRACTED_RVEF_DATA
from mecfs_bio.assets.gwas.imaging_derived_heart_phenotypes.pirruccello_et_al_2022.raw.raw_right_heart_data import \
    PIRRUCCELLO_RAW_RIGHT_HEART_DATA


def run_initial_right_heart_analysis():
    DEFAULT_RUNNER.run(
        [PIRRUCCELLO_RAW_RIGHT_HEART_DATA,
         PIRRUCCELLO_EXTRACTED_RVEF_DATA
         ],
        incremental_save=True,
    )


if __name__ == "__main__":
    run_initial_right_heart_analysis()