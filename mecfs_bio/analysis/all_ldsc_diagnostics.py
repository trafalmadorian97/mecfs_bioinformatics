"""
Materialize the LD-score-regression diagnostic plot for every GWAS on which we run an LDSC
heritability analysis.

Each plot bins a GWAS's variants by LD score and shows each bin's mean chi-square against the
fitted regression line, giving the visual diagnostic that LD-score regression normally lacks (see
ldsc_diagnostic_plot_task). LDSC_DIAGNOSTICS gathers one such task per GWAS from two sources:

- GWAS run through the standard analysis generator with a PhenotypeInfo, which now builds the
  diagnostic as a terminal task of its StandardAnalysisTaskGroup (reached here via
  ldsc_diagnostic_plot_task_unwrap, or .tasks for the add-rsids variant).
- GWAS with a standalone SNPHeritabilityByLDSCTask but no standard analysis, whose diagnostic is
  defined alongside that heritability task.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.brainstem.whole_brainstem.xue_et_al.analysis.xue_whole_brainstem_standard_analysis import (
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.c_reactive_protein.said_et_al.analysis.said_crp_standard_analysis import (
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.educational_attainment.lee_et_al_2018.analysis.lee_et_al_2018_ldsc_diagnostic_plot import (
    LEE_ET_AL_2018_LDSC_DIAGNOSTIC_PLOT,
)
from mecfs_bio.assets.gwas.height.yengo_2022.analysis.yengo_standard_analysis import (
    YENGO_HEIGHT_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.ldl.million_veterans.analysis.mv_ldl_ldsc_diagnostic_plot import (
    MV_LDL_LDSC_DIAGNOSTIC_PLOT,
)
from mecfs_bio.assets.gwas.ldl.willer_et_al.analysis.willer_ldl_standard_analysis import (
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc_diagnostic_plot import (
    DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_female.analysis.decode_me_female_standard_analysis import (
    DECODE_ME_FEMALE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_gwas_2.analysis.decode_me_gwas_2_standard_analysis import (
    DECODE_ME_GWAS_2_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_infectious_onset.analysis.decode_me_infectious_onset_standard_analysis import (
    DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_male.analysis.decode_me_male_standard_analysis import (
    DECODE_ME_MALE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_minus_pain.analysis.standard_analysis_decodeme_minus_pain_ols import (
    DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me_non_infectious_onset.analysis.decode_me_non_infectious_onset_standard_analysis import (
    DECODE_ME_NON_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.me_cfs.million_veterans.analysis.million_veterans_cfs_standard_analysis import (
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.me_cfs.neale_lab.analysis.neale_lab_cfs_standard_analysis import (
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP,
)
from mecfs_bio.assets.gwas.migraine.million_veterans.analysis.million_veterans_migraine_standard_analysis import (
    MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.migraine.uk_biobank_2025.analysis.uk_biobank_2025_migraine_standard_analysis import (
    UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_standard_analysis import (
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seronegative.analysis.ra_seronegative_standard_analysis import (
    SERONEGATIVE_RA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.analysis.ra_seropositive_standard_analysis import (
    SEROPOSITIVE_RA_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.gwas.triglycerides.willer_et_al.analysis.triglycide_standard_analysis import (
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS,
)
from mecfs_bio.build_system.task.base_task import Task

# GWAS whose diagnostic is a terminal task of a StandardAnalysisTaskGroup built directly.
_STANDARD_ANALYSIS_GROUPS = [
    XUE_WHOLE_BRAINSTEM_STANDARD_ANALYSIS,
    SAID_ET_AL_EUR_CRP_STANDARD_ANALYSIS,
    WILLER_ET_AL_EUR_LDL_STANDARD_ANALYSIS,
    WILLER_ET_AL_EUR_TG_STANDARD_ANALYSIS,
    YENGO_HEIGHT_STANDARD_ANALYSIS,
    MILLION_VETERANS_EUR_MIGRAINE_STANDARD_ANALYSIS,
    UK_BIOBANK_2025_EUR_MIGRAINE_STANDARD_ANALYSIS,
    MILLION_VETERAN_MI_EUR_STANDARD_ANALYSIS,
    MILLION_VETERANS_CFS_STANDARD_ANALYSIS_TASK_GROUP,
    NEALE_LAB_CFS_STANDARD_ANALYSIS_TASK_GROUP,
    DECODE_ME_MINUS_PAIN_OLS_STANDARD_ANALYSIS,
]

# GWAS whose standard analysis first assigns rsids; the group is wrapped one level deeper.
_STANDARD_ANALYSIS_GROUPS_WITH_RSID_ASSIGNMENT = [
    DECODE_ME_GWAS_2_STANDARD_ANALYSIS,
    DECODE_ME_MALE_STANDARD_ANALYSIS,
    DECODE_ME_FEMALE_STANDARD_ANALYSIS,
    DECODE_ME_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
    DECODE_ME_NON_INFECTIOUS_ONSET_STANDARD_ANALYSIS,
    SEROPOSITIVE_RA_STANDARD_ANALYSIS,
    SERONEGATIVE_RA_STANDARD_ANALYSIS,
]

# GWAS with a standalone LDSC heritability task and a hand-written diagnostic beside it.
_STANDALONE_DIAGNOSTICS = [
    DECODE_ME_GWAS_1_LDSC_DIAGNOSTIC_PLOT,
    MV_LDL_LDSC_DIAGNOSTIC_PLOT,
    LEE_ET_AL_2018_LDSC_DIAGNOSTIC_PLOT,
]

LDSC_DIAGNOSTICS: list[Task] = (
    [group.ldsc_diagnostic_plot_task_unwrap for group in _STANDARD_ANALYSIS_GROUPS]
    + [
        group.tasks.ldsc_diagnostic_plot_task_unwrap
        for group in _STANDARD_ANALYSIS_GROUPS_WITH_RSID_ASSIGNMENT
    ]
    + _STANDALONE_DIAGNOSTICS
)


def run_all_ldsc_diagnostics() -> None:
    """Build every LDSC diagnostic plot, saving each as it completes so a failure on one GWAS does
    not discard the plots already produced."""
    DEFAULT_RUNNER.run(LDSC_DIAGNOSTICS, incremental_save=True)


if __name__ == "__main__":
    run_all_ldsc_diagnostics()
