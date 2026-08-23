from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.kerrebeijin_ppp_rg import \
    KERREBEIJIN_ET_AL_PPP_RG_CIS_EXCLUDED
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.standard_analysis_kerrebijin_fibro import \
    KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.raw.kerrebijin_et_al_fibro_raw import \
    KERREBIJN_ET_AL_FIBRO_EUR_RAW


def go():
    DEFAULT_RUNNER.run(

            # KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS.hba_magma_tasks_unwrap.terminal_tasks()
    KERREBEIJIN_ET_AL_PPP_RG_CIS_EXCLUDED.get_terminal_tasks()

    )

if __name__ == '__main__':
    go()