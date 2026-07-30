from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_ppp_rg import (
    DECODE_ME_PPP_RG_CIS_EXCLUDED,
)


def go():
    """
    Compute genetic correlation between DecodeME and proteomic GWAS from
    the UK Biobank Pharma Proteomics Project

    """
    DEFAULT_RUNNER.run(DECODE_ME_PPP_RG_CIS_EXCLUDED.get_terminal_tasks())


if __name__ == "__main__":
    go()
