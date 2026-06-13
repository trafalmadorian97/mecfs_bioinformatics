"""
Dev script: validate the polars munge + sumstats ports against
GenomicSEM::munge / ::sumstats on the real ME/CFS + multisite-pain GWASes,
and report the speed difference.

Run with:
  pixi r python experiments/claude/runs/validate_python_munge_sumstats.py
"""

import time
from pathlib import Path

import numpy as np
import polars as pl
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genomic_sem.mecf_pain_common_factor import (
    MECFS_PAIN_COMMON_FACTOR,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    GWASMethod,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_munge import (
    munge_sumstats,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_sumstats import (
    SumstatsTrait,
    run_sumstats,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_inputs import (
    get_prevs,
    get_sample_size,
    write_munge_input,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_r_bridge import (
    run_r_munge,
)

WORKDIR = Path("experiments/claude/runs/_munge_sumstats_validation").resolve()
_METHOD_TO_FLAGS = {
    "ols": (True, False, False),
    "logistic": (False, True, False),
    "linear_prob": (False, False, True),
}


def _report(name: str, py: np.ndarray, r: np.ndarray) -> None:
    py = np.asarray(py, float)
    r = np.asarray(r, float)
    abs_err = np.max(np.abs(py - r))
    rel = np.max(np.abs(py - r) / np.maximum(np.abs(r), 1e-300))
    print(f"  {name}: max abs err {abs_err:.3e}  max rel err {rel:.3e}")


def main() -> None:
    task = MECFS_PAIN_COMMON_FACTOR
    gwas_sources = list(task.sources)
    sources = [g.source for g in gwas_sources]

    WORKDIR.mkdir(parents=True, exist_ok=True)
    tmp_dir = WORKDIR / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    munged_dir = WORKDIR / "munged_r"
    munged_dir.mkdir(parents=True, exist_ok=True)

    result = DEFAULT_RUNNER.run(
        [s.task for s in sources] + [task.ld_ref_task, task.hapmap_snps_task]
    )

    def fetch(asset_id):
        return result[asset_id]

    hapmap_path = Path(result[task.hapmap_snps_task.asset_id].path).resolve()
    sumstats_ref_path = Path(
        DEFAULT_RUNNER.run([task.sumstats_ref_task])[task.sumstats_ref_task.asset_id].path
    ).resolve()

    gsem = importr("GenomicSEM")

    input_files: list[Path] = []
    trait_names: list[str] = []
    sample_sizes: list[float] = []
    sample_prevs: list[float] = []
    pop_prevs: list[float] = []
    methods: list[GWASMethod] = []
    for g in gwas_sources:
        src = g.source
        path = write_munge_input(source=src, fetch=fetch, tmp_dir=tmp_dir)
        input_files.append(path)
        trait_names.append(src.alias)
        sample_sizes.append(get_sample_size(src))
        sp, pp = get_prevs(src.sample_info)
        sample_prevs.append(sp)
        pop_prevs.append(pp)
        methods.append(g.gwas_method)

    # ---------------- MUNGE ----------------
    print("\n==== MUNGE ====")
    t0 = time.time()
    run_r_munge(
        gsem=gsem,
        input_files=[str(p) for p in input_files],
        hapmap_path=hapmap_path,
        trait_names=trait_names,
        sample_sizes=sample_sizes,
        output_dir=munged_dir,
        info_filter=task.munge_config.info_filter,
        maf_filter=task.munge_config.maf_filter,
    )
    print(f"R munge took {time.time() - t0:.1f}s")

    ref_hm3 = pl.read_csv(hapmap_path, separator="\t")
    t0 = time.time()
    py_munged = {}
    for path, name, n in zip(input_files, trait_names, sample_sizes):
        df = pl.read_csv(path, separator="\t")
        py_munged[name] = munge_sumstats(
            df,
            ref_hm3,
            n=n,
            info_filter=task.munge_config.info_filter,
            maf_filter=task.munge_config.maf_filter,
        )
    print(f"Python munge took {time.time() - t0:.1f}s")

    for name in trait_names:
        r_m = pl.read_csv(munged_dir / f"{name}.sumstats.gz", separator="\t").sort("SNP")
        p_m = py_munged[name].sort("SNP")
        same = r_m["SNP"].to_list() == p_m["SNP"].to_list()
        print(f"[{name}] R={r_m.height} Python={p_m.height} SNPs, SNP sets equal: {same}")
        if same:
            _report(f"{name} Z", p_m["Z"].to_numpy(), r_m["Z"].to_numpy())
            print(
                f"  A1 match: {p_m['A1'].to_list() == r_m['A1'].to_list()}  "
                f"A2 match: {p_m['A2'].to_list() == r_m['A2'].to_list()}"
            )

    # ---------------- SUMSTATS ----------------
    print("\n==== SUMSTATS ====")
    se_logit = [_METHOD_TO_FLAGS[m][1] for m in methods]
    ols = [_METHOD_TO_FLAGS[m][0] for m in methods]
    linprob = [_METHOD_TO_FLAGS[m][2] for m in methods]

    t0 = time.time()
    r_ss = gsem.sumstats(
        files=ro.StrVector([str(p) for p in input_files]),
        ref=str(sumstats_ref_path),
        trait_names=ro.StrVector(trait_names),
        se_logit=ro.BoolVector(se_logit),
        OLS=ro.BoolVector(ols),
        linprob=ro.BoolVector(linprob),
        N=ro.FloatVector(sample_sizes),
        info_filter=task.sumstats_config.info_filter,
        maf_filter=task.sumstats_config.maf_filter,
        ambig=task.sumstats_config.exclude_ambig,
    )
    with localconverter(ro.default_converter + pandas2ri.converter):
        r_ss_df = pl.from_pandas(ro.conversion.get_conversion().rpy2py(r_ss))
    print(f"R sumstats took {time.time() - t0:.1f}s")

    ref_1kg = pl.read_csv(sumstats_ref_path, separator=" ")
    t0 = time.time()
    traits = [
        SumstatsTrait(
            df=pl.read_csv(path, separator="\t"),
            name=name,
            n=n,
            gwas_method=m,
        )
        for path, name, n, m in zip(input_files, trait_names, sample_sizes, methods)
    ]
    py_ss_df = run_sumstats(
        traits,
        ref_1kg,
        maf_filter=task.sumstats_config.maf_filter,
        info_filter=task.sumstats_config.info_filter,
        exclude_ambig=task.sumstats_config.exclude_ambig,
    )
    print(f"Python sumstats took {time.time() - t0:.1f}s")

    r_sorted = r_ss_df.sort("SNP")
    p_sorted = py_ss_df.sort("SNP")
    same = r_sorted["SNP"].to_list() == p_sorted["SNP"].to_list()
    print(f"R={r_sorted.height} Python={p_sorted.height} SNPs, SNP sets equal: {same}")
    if same:
        for name in trait_names:
            _report(
                f"beta.{name}",
                p_sorted[f"beta.{name}"].to_numpy(),
                r_sorted[f"beta.{name}"].to_numpy(),
            )
            _report(
                f"se.{name}",
                p_sorted[f"se.{name}"].to_numpy(),
                r_sorted[f"se.{name}"].to_numpy(),
            )


if __name__ == "__main__":
    main()
