"""
Analyze diastolic blood pressure  GWAS
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.analysis.keaton_dbp_standard_analysis import \
    KEATON_DBP_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.blood_pressure.keaton_et_al_diastolic.raw.raw_dbp_data import KEATON_ET_AL_DBP_RAW
from mecfs_bio.assets.gwas.red_blood_cell_volume.million_veterans.analysis.rbc_vol_sldsc import (
    MILLION_VETERANS_EUR_RBC_VOL_S_LDSC_TASKS,
)
from mecfs_bio.assets.gwas.red_blood_cell_volume.million_veterans.analysis.rbc_volume_magma_task_generator import (
    MILLION_VETERANS_EUR_RBC_VOLUME_MAGMA_TASKS,
)


def run_dbp_analysis():
    """
    Analyze Keaton meta GWAS f diastolic blood pressure

    """
    DEFAULT_RUNNER.run(
        [
            KEATON_ET_AL_DBP_RAW
        ]+KEATON_DBP_STANDARD_ANALYSIS.get_terminal_tasks() ,
        incremental_save=True,
        )


if __name__ == "__main__":
    run_dbp_analysis()
