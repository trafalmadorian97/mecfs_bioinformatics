"""
Sweep the synthetic problem size of the SuSiE unit test fixture.

For each (n, m) configuration and each random seed, simulate data exactly as the
susie_prerequisite_file_tasks fixture in
test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py does, run the real
SusieRFinemapTask, and record:

- wall time for data simulation and for the task run,
- whether the test_fine_mapping assertions would pass (PIP >= 0.95 at every causal
  variant, and number of credible sets == number of causal variants),
- the margin: min PIP over causal variants (causal case), or the number of credible
  sets and max PIP (null case).

n = number of simulated individuals; m = number of variants in the locus.
Causal configurations: "two" puts unit effects on the first and last variant, "null"
has no causal variant (both match the fixture params).

Usage: pixi r python experiments/claude/susie_test_speed/susie_problem_size_sweep.py \
    | tee experiments/claude/susie_test_speed/sweep.log
"""

import sys
import tempfile
import time
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
import rpy2.robjects as ro
import scipy.sparse
from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.simple_meta_to_path import (
    SimpleMetaToPath,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.simple_hasher import (
    SimpleHasher,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_rebuilder_core import (
    VerifyingTraceRebuilder,
)
from mecfs_bio.build_system.scheduler.topological_scheduler import topological
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    CS_DATA_SUBDIR,
    PIP_COLUMN,
    PIP_FILENAME,
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
)
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.util.r_packages.cached_importr import cached_importr

CONFIGS = [(2500, 100), (1000, 100), (500, 100), (2500, 50), (1000, 50), (500, 50), (500, 30)]
N_SEEDS = 25
FIXTURE_SEED = 40


def simulate_and_write(n: int, m: int, causal: list[int], seed: int, scratch: Path) -> None:
    """Fixture data generation, copied from susie_prerequisite_file_tasks."""
    generator = np.random.default_rng(seed)
    true_effects = np.zeros(m)
    for cv in causal:
        true_effects[cv] = 1
    q = np.linspace(1, 0, num=m)
    lamb = 0.2
    covar = (1 - lamb) * q.reshape(-1, 1) * q.reshape(1, -1) + lamb * np.eye(m)
    genotypes = generator.multivariate_normal(np.zeros(m), cov=covar, size=n)
    phenotypes = genotypes @ true_effects + generator.normal(loc=0, scale=0.1, size=n)
    conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
    gwas = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1] * m,
            GWASLAB_POS_COL: list(range(m)),
            GWASLAB_EFFECT_ALLELE_COL: "A",
            GWASLAB_NON_EFFECT_ALLELE_COL: "C",
        }
    )
    ref = gwas.copy().iloc[::-1]
    ld = np.corrcoef(genotypes.T)[::-1, ::-1]
    with localconverter(conv):
        g_r = ro.conversion.get_conversion().py2rpy(genotypes)
        p_r = ro.conversion.get_conversion().py2rpy(phenotypes)
    result = cached_importr("susieR").univariate_regression(g_r, p_r)
    with localconverter(conv):
        gwas[GWASLAB_BETA_COL] = result.rx2("betahat")
        gwas[GWASLAB_SE_COL] = result.rx2("sebetahat")
    gwas.to_parquet(scratch / "gwas_data")
    ref.to_parquet(scratch / "ld_labels")
    scipy.sparse.save_npz(scratch / "ld_matrix.npz", scipy.sparse.coo_matrix(ld / 2))


def run_task(n: int, scratch: Path) -> Path:
    parquet = DataFrameReadSpec(DataFrameParquetFormat())
    susie = SusieRFinemapTask(
        meta=SimpleDirectoryMeta(AssetId("directory")),
        gwas_data_task=ExternalFileCopyTask(
            SimpleFileMeta(AssetId("gwas_data"), read_spec=parquet),
            external_path=scratch / "gwas_data",
        ),
        ld_labels_task=ExternalFileCopyTask(
            SimpleFileMeta("ld_labels", read_spec=parquet),
            external_path=scratch / "ld_labels",
        ),
        ld_matrix_source=BroadInstituteFormatLDMatrix(
            ExternalFileCopyTask(
                SimpleFileMeta("ld_matrix"), external_path=scratch / "ld_matrix.npz"
            )
        ),
        effective_sample_size=n,
        max_credible_sets=10,
    )
    asset_dir = scratch / "asset_dir"
    asset_dir.mkdir()
    store, _ = topological(
        rebuilder=VerifyingTraceRebuilder(SimpleHasher.md5_hasher()),
        tasks=find_tasks([susie]),
        targets=[susie.asset_id],
        wf=make_wf(),
        info=VerifyingTraceInfo.empty(),
        meta_to_path=SimpleMetaToPath(root=asset_dir),
    )
    asset = store[susie.asset_id]
    assert isinstance(asset, DirectoryAsset)
    return asset.path


def evaluate(n: int, m: int, causal_name: str, seed: int) -> dict[str, object]:
    causal = [0, m - 1] if causal_name == "two" else []
    with tempfile.TemporaryDirectory() as d:
        scratch = Path(d)
        t0 = time.perf_counter()
        simulate_and_write(n=n, m=m, causal=causal, seed=seed, scratch=scratch)
        t1 = time.perf_counter()
        out = run_task(n=n, scratch=scratch)
        t2 = time.perf_counter()
        pip = pd.read_parquet(out / PIP_FILENAME)[PIP_COLUMN].to_numpy()
        n_cs = len(list((out / CS_DATA_SUBDIR).glob("*")))
    min_causal_pip = float(min(pip[c] for c in causal)) if causal else float("nan")
    passed = n_cs == len(causal) and all(pip[c] >= 0.95 for c in causal)
    return {
        "n": n,
        "m": m,
        "causal": causal_name,
        "seed": seed,
        "sim_s": t1 - t0,
        "task_s": t2 - t1,
        "n_cs": n_cs,
        "min_causal_pip": min_causal_pip,
        "max_pip": float(pip.max()),
        "passed": passed,
    }


def main() -> None:
    evaluate(2500, 100, "two", FIXTURE_SEED)  # warm up R package loading
    seeds = [FIXTURE_SEED] + [1000 + i for i in range(N_SEEDS - 1)]
    rows = []
    for n, m in CONFIGS:
        for causal_name in ["two", "null"]:
            for seed in seeds:
                rows.append(evaluate(n, m, causal_name, seed))
            print(f"done n={n} m={m} causal={causal_name}", file=sys.stderr, flush=True)
    df = pl.DataFrame(rows)
    df.write_csv(Path(__file__).parent / "sweep_results.csv")
    summary = (
        df.group_by("n", "m", "causal", maintain_order=True)
        .agg(
            pl.col("passed").mean().alias("pass_rate"),
            pl.col("passed").filter(pl.col("seed") == FIXTURE_SEED).first().alias("fixture_seed_passes"),
            pl.col("min_causal_pip").min().alias("worst_min_causal_pip"),
            pl.col("max_pip").max().alias("worst_null_max_pip"),
            pl.col("sim_s").median().alias("median_sim_s"),
            pl.col("task_s").median().alias("median_task_s"),
        )
    )
    with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_width_chars=200):
        print(summary)


if __name__ == "__main__":
    main()
