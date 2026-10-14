"""Build the hg19 genetic-map parquet and report its shape."""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genetic_map.genetic_map_hg19 import (
    GENETIC_MAP_HG19,
)


def main() -> None:
    store = DEFAULT_RUNNER.run([GENETIC_MAP_HG19], incremental_save=True)
    p = store[GENETIC_MAP_HG19.asset_id].path
    df = pl.read_parquet(p)
    print("OUTPUT:", p)
    print("rows:", df.height, "cols:", df.columns)
    print("chroms:", sorted(df["CHR"].unique().to_list()))
    print(df.head(3))
    print("rate stats:", df["recomb_rate_cm_per_mb"].describe())


if __name__ == "__main__":
    main()
