"""
Materialize the GRCh38 genome-reference harmonized trait tasks added for the CSF coloc work.

Usage: pixi r python experiments/claude/csf_coloc/materialize_traits.py <ad|pd|scz|als>...
"""

import sys

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.als.van_rheenen_et_al_2021.processed.van_rheenen_2021_build_38_harmonized import (
    VAN_RHEENEN_2021_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.processed.bellenguez_build_38_harmonized import (
    BELLENGUEZ_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.parkinsons.nalls_et_al_2019.processed.nalls_2019_build_38_harmonized import (
    NALLS_2019_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.schizophrenia.pgc2022.processed.pgc_2022_sch_build_38_harmonized import (
    PGC_2022_SCH_BUILD_38_HARMONIZED,
)
from mecfs_bio.build_system.task.base_task import Task

TASKS: dict[str, Task] = {
    "ad": BELLENGUEZ_BUILD_38_HARMONIZED,
    "pd": NALLS_2019_BUILD_38_HARMONIZED,
    "scz": PGC_2022_SCH_BUILD_38_HARMONIZED,
    "als": VAN_RHEENEN_2021_BUILD_38_HARMONIZED,
}


def main() -> None:
    for name in sys.argv[1:]:
        DEFAULT_RUNNER.run(
            [TASKS[name]], incremental_save=True, must_rebuild_transitive=[]
        )


if __name__ == "__main__":
    main()
