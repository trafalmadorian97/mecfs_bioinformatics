from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.bmi.pulit_et_al_2018.analysis.pulit_bmi_ppp_rg import (
    PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED,
)


def go() -> None:
    """
    Run genetic correlation analysis via CT-LDSC correlating the GIANT + UK Biobank GWAS of body mass index
    (Pulit et al. 2019) and plasma protein levels from the UK Biobank Pharma Proteomics Project.
    """
    DEFAULT_RUNNER.run(PULIT_ET_AL_2018_BMI_PPP_RG_CIS_EXCLUDED.get_terminal_tasks())


if __name__ == "__main__":
    go()
