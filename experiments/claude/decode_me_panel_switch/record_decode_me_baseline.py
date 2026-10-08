"""
Snapshot the key built DecodeME outputs that the switch to the gnomAD panel will rebuild.

Commit 76d81a00 switched rsID-assignment harmonization from 1000 Genomes EUR to the gnomAD
v2.1.1 nfe_nwe panel (and eb8c7a31 added the palindrome distance check). The next build will
rebuild every asset downstream of a DecodeME genome-reference harmonization. This script must
run BEFORE that build: it never runs the build system, it only reads the asset store.

Scope: the analyses in ROOT_MODULES and every module in ROOT_PACKAGES (fine mapping). Within
the task graphs of those analyses, the recorded tasks are the DecodeME
GenomeReferenceHarmonizationTasks and every task downstream of one; tasks upstream of
harmonization (raw downloads, reference data) cannot change and are skipped.

For each recorded task whose output exists on disk: its path, size and per-file md5 go to
manifest.json, and the output is copied into the baseline directory when it is at most
MAX_COPY_BYTES. Larger outputs (the harmonized and rsID-joined tables) are hashed only.

The baseline lives under assets/ (gitignored), so it is never committed.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.record_decode_me_baseline \
        2>&1 | tee experiments/claude/decode_me_panel_switch/record_decode_me_baseline.log
"""

import hashlib
import importlib
import json
import pkgutil
import shutil
from collections.abc import Iterator, Sequence
from pathlib import Path
from types import ModuleType

from attrs import asdict, frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks

_ANALYSIS = "mecfs_bio.assets.gwas.me_cfs.decode_me.analysis"
ROOT_MODULES = (
    f"{_ANALYSIS}.decode_me_gwas_1_cis_pqtl_mr",
    f"{_ANALYSIS}.decode_me_gwas_1_ldsc",
    f"{_ANALYSIS}.decode_me_sldsc",
    f"{_ANALYSIS}.magma.decode_me_hba_magma_analysis",
    f"{_ANALYSIS}.magma.decode_me_gwas_1_build_37_magma_ensembl_specific_tissue_gene_sets",
)
ROOT_PACKAGES = (
    f"{_ANALYSIS}.fine_mapping.with_palindromes",
    f"{_ANALYSIS}.fine_mapping.polyfun_explainability",
    f"{_ANALYSIS}.fine_mapping.polyfun_approach_2_l2_sldsc_prior",
)
BASELINE_DIR = Path("assets") / "baselines" / "decode_me_gnomad_panel_switch"
MANIFEST_NAME = "manifest.json"
AFFECTED_HARMONIZATION_PREFIX = "decode_me"
MAX_COPY_BYTES = 500 * 1024**2
_TERMINAL_METHODS = ("terminal_tasks", "get_terminal_tasks")
_HASH_CHUNK_BYTES = 16 * 1024**2


@frozen(slots=True)
class AssetRecord:
    asset_id: str
    task_type: str
    path: str
    exists: bool
    size_bytes: int
    copied: bool
    file_md5: dict[str, str]


def module_tasks(value: object) -> Iterator[Task]:
    """Tasks held by a module attribute: a Task, a task group, or a list or dict of them."""
    if isinstance(value, type):
        return
    if isinstance(value, Task):
        yield value
        return
    if isinstance(value, (list, tuple)):
        for item in value:
            yield from module_tasks(item)
        return
    if isinstance(value, dict):
        for item in value.values():
            yield from module_tasks(item)
        return
    for method in _TERMINAL_METHODS:
        terminal = getattr(value, method, None)
        if callable(terminal):
            yield from terminal()
            return


def root_modules() -> list[ModuleType]:
    modules = [importlib.import_module(name) for name in ROOT_MODULES]
    for package_name in ROOT_PACKAGES:
        package = importlib.import_module(package_name)
        modules += [
            importlib.import_module(info.name)
            for info in pkgutil.iter_modules(package.__path__, prefix=package_name + ".")
        ]
    return modules


def root_tasks(modules: Sequence[ModuleType]) -> list[Task]:
    """One task per asset id: modules re-export the same tasks, and find_tasks compares
    duplicates by equality, which fails for tasks holding polars expressions."""
    unique: dict[AssetId, Task] = {}
    for module in modules:
        for value in vars(module).values():
            for task in module_tasks(value):
                unique.setdefault(task.asset_id, task)
    return list(unique.values())


def affected_task_ids(tasks: dict[AssetId, Task]) -> set[AssetId]:
    """DecodeME harmonizations and every task that depends on one, directly or not."""
    affected = {
        asset_id
        for asset_id, task in tasks.items()
        if isinstance(task, GenomeReferenceHarmonizationTask)
        and str(asset_id).startswith(AFFECTED_HARMONIZATION_PREFIX)
    }
    assert affected, "no DecodeME genome-reference harmonization tasks found"
    changed = True
    while changed:
        changed = False
        for asset_id, task in tasks.items():
            if asset_id not in affected and any(
                dep.asset_id in affected for dep in task.deps
            ):
                affected.add(asset_id)
                changed = True
    return affected


def files_under(path: Path) -> list[Path]:
    if path.is_file():
        return [path]
    return sorted(p for p in path.rglob("*") if p.is_file())


def md5_of(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(_HASH_CHUNK_BYTES), b""):
            digest.update(chunk)
    return digest.hexdigest()


def record_asset(task: Task, path: Path) -> AssetRecord:
    if not path.exists():
        return AssetRecord(
            asset_id=str(task.asset_id),
            task_type=type(task).__name__,
            path=str(path),
            exists=False,
            size_bytes=0,
            copied=False,
            file_md5={},
        )
    files = files_under(path)
    size = sum(f.stat().st_size for f in files)
    root = path if path.is_dir() else path.parent
    copied = size <= MAX_COPY_BYTES
    if copied:
        target = BASELINE_DIR / "outputs" / str(task.asset_id)
        for f in files:
            destination = target / f.relative_to(root)
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(f, destination)
    return AssetRecord(
        asset_id=str(task.asset_id),
        task_type=type(task).__name__,
        path=str(path),
        exists=True,
        size_bytes=size,
        copied=copied,
        file_md5={str(f.relative_to(root)): md5_of(f) for f in files},
    )


def report(records: Sequence[AssetRecord], n_tasks: int) -> None:
    built = [r for r in records if r.exists]
    print(f"tasks in the graphs of the chosen analyses: {n_tasks}")
    print(f"downstream of a DecodeME harmonization: {len(records)}")
    print(f"  built and recorded: {len(built)}")
    print(f"  copied: {sum(r.copied for r in built)}")
    print(f"  hashed only (over {MAX_COPY_BYTES // 1024**2} MiB): {sum(not r.copied for r in built)}")
    print(f"  not built (nothing to record): {len(records) - len(built)}")
    print(f"bytes recorded: {sum(r.size_bytes for r in built):,}")
    for r in records:
        status = "copied" if r.copied else ("hashed" if r.exists else "NOT BUILT")
        print(f"  {status:9} {r.task_type:40} {r.asset_id}")


def main() -> None:
    manifest = BASELINE_DIR / MANIFEST_NAME
    assert not manifest.exists(), (
        f"a baseline already exists at {BASELINE_DIR}; move it aside before re-recording"
    )
    tasks = dict(find_tasks(root_tasks(root_modules())).tasks)
    meta_to_path = DEFAULT_RUNNER.meta_to_path
    records = [
        record_asset(tasks[asset_id], meta_to_path(tasks[asset_id].meta))
        for asset_id in sorted(affected_task_ids(tasks), key=str)
    ]
    BASELINE_DIR.mkdir(parents=True, exist_ok=True)
    manifest.write_text(json.dumps([asdict(r) for r in records], indent=2) + "\n")
    report(records, n_tasks=len(tasks))
    print(f"manifest: {manifest}")


if __name__ == "__main__":
    main()
