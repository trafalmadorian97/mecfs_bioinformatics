"""
Real-data check of the fast all-pairs CT-LDSC generator on the initial study traits.

Two modes:

- diagnostics: build only the aligned-trait tasks, with permissive guards, and tabulate each
  trait's alignment diagnostics (coverage, allele_mismatch_fraction, position_concordance) and
  where its declared genome build came from (source meta, the gwaslab pickle, or nowhere). These
  numbers set the AlignmentGuards defaults.
- compare: run the fast generator end to end (timed), build the existing gwaslab generator's
  pairs and heritability aggregations (slow when not already cached), and tabulate per pair
  rg_fast, rg_gwaslab, delta rg and delta rg / se_gwaslab, per trait h2 and delta h2, and which
  pairs pass Bonferroni under each method.

Both modes write to the shared asset store configured in default_runner_config.yaml.

Run: pixi r python experiments/claude/fast_ct_ldsc_vs_gwaslab/compare.py diagnostics|compare
Output is tee'd to compare_<mode>.log next to this script.
"""

from __future__ import annotations

import sys
import time
from pathlib import Path
from typing import Literal

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.fast_genetic_correlation_asset_generator import (
    fast_genetic_corr_by_ct_ldsc_asset_generator,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_fast_initial_asset_generator import (
    CT_LDSC_FAST_INITIAL_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_asset_generator import (
    CT_LDSC_INITIAL_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_sources import (
    CT_LDSC_INITIAL_SOURCES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fast_ct_ldsc.aligned_trait_task import (
    read_aligned_trait,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import AlignmentGuards

Mode = Literal["diagnostics", "compare"]
_HERE = Path(__file__).resolve().parent
_PERMISSIVE = AlignmentGuards(
    min_coverage=0.0, max_allele_mismatch_fraction=1.0, min_position_concordance=0.0
)
_ALPHA = 0.05


class _Tee:
    def __init__(self, path: Path) -> None:
        self._file = path.open("w")
        self._stdout = sys.stdout

    def write(self, text: str) -> int:
        self._file.write(text)
        self._file.flush()
        return self._stdout.write(text)

    def flush(self) -> None:
        self._file.flush()
        self._stdout.flush()


def _fetcher(built: dict[AssetId, Asset]):
    def fetch(asset_id: AssetId) -> Asset:
        return built[asset_id]

    return fetch


def _read_table(task: Task, built: dict[AssetId, Asset]) -> pl.DataFrame:
    asset = built[task.asset_id]
    assert isinstance(asset, FileAsset)
    return (
        scan_dataframe_asset(asset, meta=task.meta, parquet_backend="polars")
        .collect(backend="polars")
        .to_polars()
    )


def diagnostics() -> None:
    gen = fast_genetic_corr_by_ct_ldsc_asset_generator(
        "initial_rg_fast_diagnostics", sources=CT_LDSC_INITIAL_SOURCES, guards=_PERMISSIVE
    )
    tasks = list(gen.aligned_trait_tasks.values())
    start = time.perf_counter()
    built = dict(DEFAULT_RUNNER.run(tasks, incremental_save=True))
    print(f"aligned {len(tasks)} traits in {time.perf_counter() - start:.1f} s")
    rows = []
    for source, task in zip(
        sorted(CT_LDSC_INITIAL_SOURCES, key=lambda s: s.alias.lower()), tasks
    ):
        diag = read_aligned_trait(task, _fetcher(built)).diagnostics
        meta_build = getattr(source.task.meta, "genome_build", None)
        build_source = (
            "meta"
            if meta_build is not None
            else ("pickle" if diag.declared_build is not None else "none")
        )
        rows.append(
            {
                "alias": task.alias,
                "n_context": diag.n_context,
                "n_found": diag.n_found,
                "n_matched": diag.n_matched,
                "coverage": diag.coverage,
                "allele_mismatch_fraction": diag.allele_mismatch_fraction,
                "position_concordance": diag.position_concordance,
                "declared_build": diag.declared_build,
                "build_source": build_source,
                "source_meta": type(source.task.meta).__name__,
            }
        )
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=250):
        print(pl.DataFrame(rows))


def compare() -> None:
    fast = CT_LDSC_FAST_INITIAL_ASSET_GENERATOR
    start = time.perf_counter()
    built = dict(DEFAULT_RUNNER.run(fast.terminal_tasks(), incremental_save=True))
    print(f"fast generator end to end: {time.perf_counter() - start:.1f} s")
    old = CT_LDSC_INITIAL_ASSET_GENERATOR
    start = time.perf_counter()
    built.update(
        DEFAULT_RUNNER.run(
            [old.aggregation_task, old.heritability_aggregation_task],
            incremental_save=True,
        )
    )
    print(f"gwaslab aggregations (cached or built): {time.perf_counter() - start:.1f} s")

    fast_pairs = _read_table(fast.pairs_task, built)
    old_pairs = _read_table(old.aggregation_task, built).select(
        pl.col("p1").cast(pl.String),
        pl.col("p2").cast(pl.String),
        pl.col("rg").cast(pl.Float64).alias("rg_gwaslab"),
        pl.col("se").cast(pl.Float64).alias("se_gwaslab"),
        pl.col("p").cast(pl.Float64).alias("p_gwaslab"),
    )
    n_pairs = fast_pairs.height
    threshold = _ALPHA / n_pairs
    pairs = (
        fast_pairs.select(
            "p1",
            "p2",
            pl.col("rg").alias("rg_fast"),
            pl.col("se").alias("se_fast"),
            pl.col("p").alias("p_fast"),
        )
        .join(old_pairs, on=["p1", "p2"], how="full", coalesce=True)
        .with_columns(
            (pl.col("rg_fast") - pl.col("rg_gwaslab")).alias("d_rg"),
        )
        .with_columns(
            (pl.col("d_rg") / pl.col("se_gwaslab")).alias("d_rg_over_se"),
            (pl.col("p_fast") < threshold).alias("bonf_fast"),
            (pl.col("p_gwaslab") < threshold).alias("bonf_gwaslab"),
        )
        .sort(pl.col("d_rg_over_se").abs(), descending=True, nulls_last=True)
    )
    fast_h2 = _read_table(fast.heritability_task, built)
    old_h2 = _read_table(old.heritability_aggregation_task, built)
    h2 = fast_h2.select(
        "p",
        pl.col("h2_obs").alias("h2_obs_fast"),
        pl.col("h2_obs_se").alias("h2_obs_se_fast"),
        pl.col("Intercept").alias("intercept_fast"),
    ).join(
        old_h2.select(
            pl.col("p").cast(pl.String),
            pl.col("h2_obs").cast(pl.Float64).alias("h2_obs_gwaslab"),
            pl.col("h2_obs_se").cast(pl.Float64).alias("h2_obs_se_gwaslab"),
            pl.col("Intercept").cast(pl.Float64).alias("intercept_gwaslab"),
        ),
        on="p",
        how="full",
        coalesce=True,
    ).with_columns(
        (pl.col("h2_obs_fast") - pl.col("h2_obs_gwaslab")).alias("d_h2"),
        (
            (pl.col("h2_obs_fast") - pl.col("h2_obs_gwaslab"))
            / pl.col("h2_obs_se_gwaslab")
        ).alias("d_h2_over_se"),
    )
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=250, float_precision=4):
        print(pairs)
        print(h2)
    changed = pairs.filter(pl.col("bonf_fast") != pl.col("bonf_gwaslab"))
    print(f"Bonferroni threshold {threshold:.3g} over {n_pairs} pairs")
    print(f"max |d_rg / se_gwaslab| = {pairs['d_rg_over_se'].abs().max()}")
    print(f"pairs whose Bonferroni call differs: {changed.height}")
    print(changed.select("p1", "p2", "p_fast", "p_gwaslab"))


def main(mode: Mode) -> None:
    sys.stdout = _Tee(_HERE / f"compare_{mode}.log")
    if mode == "diagnostics":
        diagnostics()
    else:
        compare()


if __name__ == "__main__":
    arg = sys.argv[1] if len(sys.argv) > 1 else ""
    assert arg in ("diagnostics", "compare"), "usage: compare.py diagnostics|compare"
    main("diagnostics" if arg == "diagnostics" else "compare")
