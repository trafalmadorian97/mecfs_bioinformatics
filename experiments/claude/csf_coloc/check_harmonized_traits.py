"""
Sanity checks of the GRCh38 harmonized trait tasks added for the CSF coloc work: row count,
N range, genome-wide-significant variant count, the top variants in named regions, and the
position of named rsIDs (checked by hand against Ensembl GRCh38).

Usage: pixi r python experiments/claude/csf_coloc/check_harmonized_traits.py <ad|pd|scz|als>...
"""

import sys

import polars as pl

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
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

# (task, regions (chrom, start, end) expected to hold GWS signal, rsIDs whose position to print)
CHECKS: dict[str, tuple[Task, list[tuple[int, int, int]], list[str]]] = {
    "ad": (
        BELLENGUEZ_BUILD_38_HARMONIZED,
        [(19, 44_400_000, 45_400_000)],
        ["rs429358"],
    ),
    "pd": (NALLS_2019_BUILD_38_HARMONIZED, [(4, 89_200_000, 90_200_000)], ["rs356182"]),
    "scz": (
        PGC_2022_SCH_BUILD_38_HARMONIZED,
        [(6, 25_000_000, 34_000_000)],
        ["rs13107325"],
    ),
    "als": (
        VAN_RHEENEN_2021_BUILD_38_HARMONIZED,
        [(9, 27_000_000, 28_000_000), (19, 17_200_000, 18_000_000)],
        ["rs3849943", "rs12608932"],
    ),
}
_Z = (pl.col(GWASLAB_BETA_COL) / pl.col(GWASLAB_SE_COL)).abs()
_GWS_Z = 5.45


def check(name: str) -> None:
    task, regions, rsids = CHECKS[name]
    frame = pl.scan_parquet(DEFAULT_RUNNER.meta_to_path(task.meta))
    cols = frame.collect_schema().names()
    summary = frame.select(
        pl.len().alias("rows"),
        pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).min().alias("n_min"),
        pl.col(GWASLAB_SAMPLE_SIZE_COLUMN).max().alias("n_max"),
        (_Z > _GWS_Z).sum().alias("gws_variants"),
    ).collect()
    print(f"== {name} {task.asset_id}\n{summary}")
    for chrom, start, end in regions:
        top = (
            frame.filter(
                pl.col(GWASLAB_CHROM_COL) == chrom,
                pl.col(GWASLAB_POS_COL).is_between(start, end),
            )
            .with_columns(_Z.alias("absZ"))
            .sort("absZ", descending=True)
            .head(3)
            .collect()
        )
        print(f"top in chr{chrom}:{start}-{end}\n{top}")
    if GWASLAB_RSID_COL in cols:
        print(frame.filter(pl.col(GWASLAB_RSID_COL).is_in(rsids)).collect())
    else:
        print(f"no {GWASLAB_RSID_COL} column; rsID positions not checked here")


def main() -> None:
    for name in sys.argv[1:]:
        check(name)


if __name__ == "__main__":
    main()
