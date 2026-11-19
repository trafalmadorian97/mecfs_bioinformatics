from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.analysis.kerrebijn_fibro_ppp_coloc import (
    KERREBIJN_FIBRO_PPP_COLOC_FULL,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_pain_ppp_coloc import (
    JOHNSTON_PAIN_PPP_COLOC_FULL,
)


def go():
    DEFAULT_RUNNER.run(
        [
            KERREBIJN_FIBRO_PPP_COLOC_FULL.display_task,
            JOHNSTON_PAIN_PPP_COLOC_FULL.display_task,
        ]
    )


if __name__ == "__main__":
    go()
