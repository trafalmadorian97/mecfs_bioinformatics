"""
Dev script: validate the pure-Python LDSC (`_genomic_sem_ldsc.run_ldsc`)
against GenomicSEM::ldsc on the real ME/CFS + multisite-pain GWASes.

It materializes the dependencies of the common-factor example, munges the two
traits in R (shared step), runs both R `ldsc` and the Python port on the same
munged files + LD reference, and prints the element-wise agreement of S/V/I.

Run with:  pixi r python experiments/tralfamadorian97/runs/validate_python_ldsc.py
"""

from pathlib import Path

import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import importr

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.multi_trait.genomic_sem.mecf_pain_common_factor import (
    MECFS_PAIN_COMMON_FACTOR,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_ldsc import run_ldsc
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    _get_prevs,
    _get_sample_size,
    _ld_dir_with_genomic_sem_naming,
    _run_ldsc,
    _run_munge,
    _write_munge_input,
)

WORKDIR = Path("experiments/tralfamadorian97/runs/_ldsc_validation").resolve()
LD_PREFIX = "LDscore."


def _r_mat(covstruc, name: str) -> np.ndarray:
    conv = ro.default_converter + pandas2ri.converter
    with localconverter(conv):
        return np.asarray(ro.conversion.get_conversion().rpy2py(covstruc.rx2(name)))


def main() -> None:
    task = MECFS_PAIN_COMMON_FACTOR
    gwas_sources = list(task.sources)
    sources = [g.source for g in gwas_sources]

    WORKDIR.mkdir(parents=True, exist_ok=True)
    munged_dir = WORKDIR / "munged"
    munged_dir.mkdir(parents=True, exist_ok=True)
    tmp_dir = WORKDIR / "tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # 1. Materialize raw GWAS inputs + reference assets.
    dep_tasks = (
        [s.task for s in sources] + [task.ld_ref_task, task.hapmap_snps_task]
    )
    print("Materializing dependencies ...")
    result = DEFAULT_RUNNER.run(dep_tasks)

    def fetch(asset_id):
        return result[asset_id]

    # munge chdirs into the output dir, so all paths handed to R must be absolute.
    ld_path = str(Path(result[task.ld_ref_task.asset_id].path).resolve())
    hapmap_path = str(Path(result[task.hapmap_snps_task.asset_id].path).resolve())

    gsem = importr("GenomicSEM")

    # 2. Write munge inputs and run R munge (shared step).
    input_files: list[str] = []
    trait_names: list[str] = []
    sample_sizes: list[float] = []
    sample_prevs: list[float] = []
    population_prevs: list[float] = []
    for g in gwas_sources:
        src = g.source
        path = _write_munge_input(source=src, fetch=fetch, tmp_dir=tmp_dir)
        input_files.append(str(path))
        trait_names.append(src.alias)
        sample_sizes.append(_get_sample_size(src))
        sp, pp = _get_prevs(src.sample_info)
        sample_prevs.append(sp)
        population_prevs.append(pp)

    print("Running R munge ...")
    _run_munge(
        gsem=gsem,
        input_files=input_files,
        hapmap_path=hapmap_path,
        trait_names=trait_names,
        sample_sizes=sample_sizes,
        output_dir=munged_dir,
        info_filter=task.munge_config.info_filter,
        maf_filter=task.munge_config.maf_filter,
    )
    munged_paths = [munged_dir / f"{name}.sumstats.gz" for name in trait_names]
    for p in munged_paths:
        assert p.is_file(), f"missing munged file {p}"

    # 3. R ldsc (ground truth).
    print("Running R ldsc ...")
    covstruc = _run_ldsc(
        gsem=gsem,
        munged_paths=[str(p) for p in munged_paths],
        trait_names=trait_names,
        sample_prevs=sample_prevs,
        population_prevs=population_prevs,
        ld_path=ld_path,
        ld_file_basename_prefix=LD_PREFIX,
        ldsc_log_prefix=str(WORKDIR / "r_ldsc"),
    )
    r_S = _r_mat(covstruc, "S")
    r_V = _r_mat(covstruc, "V")
    r_I = _r_mat(covstruc, "I")
    r_S_stand = _r_mat(covstruc, "S_Stand")

    # 4. Python ldsc on the same munged files + LD reference.
    print("Running Python run_ldsc ...")
    with _ld_dir_with_genomic_sem_naming(ld_path, LD_PREFIX) as eff_ld:
        py = run_ldsc(
            munged_paths=munged_paths,
            ld_dir=Path(eff_ld),
            sample_prev=sample_prevs,
            population_prev=population_prevs,
            n_chrom=22,
            stand=True,
        )

    # 5. Report.
    def report(name: str, r: np.ndarray, p: np.ndarray) -> None:
        r = np.asarray(r, dtype=float)
        p = np.asarray(p, dtype=float)
        abs_err = np.max(np.abs(r - p))
        denom = np.maximum(np.abs(r), 1e-300)
        rel_err = np.max(np.abs(r - p) / denom)
        print(f"\n=== {name} ===")
        print("R:\n", r)
        print("Python:\n", p)
        print(f"max abs err = {abs_err:.3e}   max rel err = {rel_err:.3e}")

    print("\n================ COMPARISON ================")
    print("trait order:", trait_names)
    report("S", r_S, py.S)
    report("I", r_I, py.I)
    report("V", r_V, py.V)
    report("S_Stand", r_S_stand, py.S_Stand)
    print(f"\nN (R):      {np.asarray(covstruc.rx2('N')).ravel()}")
    print(f"N (Python): {py.N}")
    print(f"m (R): {float(np.asarray(covstruc.rx2('m')))}   m (Python): {py.m}")


if __name__ == "__main__":
    main()
