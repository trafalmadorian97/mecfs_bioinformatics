from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_pp_rg import BELLENGUEZ_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.raw.raw_bellenguez_data import BELLENGUEZ_ET_AL_ALZHEIMERS_RAW


def go():
    DEFAULT_RUNNER.run(
        BELLENGUEZ_PPP_RG_CIS_EXCLUDED.get_terminal_tasks()

    )

if __name__ == '__main__':
    go()