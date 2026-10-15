"""
Script to analyze the GIANT + UK Biobank GWAS of body mass index from Pulit et al.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_standard_analysis import (
    PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS,
)


def run_bmi_analysis() -> None:
    """
    Apply the standard analysis pipeline to the BMI GWAS of Pulit et al. (2019).
    Analysis includes:
    - Application of MAGMA to GTEx data to identify possible key tissues.
    - Application of MAGMA to human brain atlas to identify key brain tissue
    - Use of S-LDSC to investigate key cell and tissue types
    - Estimation of SNP heritability via LDSC
    """
    DEFAULT_RUNNER.run(
        PULIT_ET_AL_2018_BMI_STANDARD_ANALYSIS.get_terminal_tasks(),
        incremental_save=True,
        must_rebuild_transitive=[],
    )


if __name__ == "__main__":
    run_bmi_analysis()
