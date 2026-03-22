from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.human_herpesvirus_7_dna.kamitaki_et_al_2025.analysis.kamitaki_et_al_2025_standard_analysis import \
    KAMITAKI_ET_AL_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.human_herpesvirus_7_dna.kamitaki_et_al_2025.raw.kamitaki_et_al_2025_raw import \
    KAMITAKI_ET_AL_HHV7_DNA_RAW


def run_hhv7():
    DEFAULT_RUNNER.run(
        [
            KAMITAKI_ET_AL_HHV7_DNA_RAW,
        ] + KAMITAKI_ET_AL_STANDARD_ANALYSIS.terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[
        ]
    )


if __name__ == "__main__":
    run_hhv7()
