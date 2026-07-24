from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.myocardial_infarction.analysis.mi_ppp_rg import MV_MI_PPP_RG_CIS_EXCLUDED


def go():
    DEFAULT_RUNNER.run(

        MV_MI_PPP_RG_CIS_EXCLUDED.get_terminal_tasks()

    )

if __name__ == "__main__":
    go()