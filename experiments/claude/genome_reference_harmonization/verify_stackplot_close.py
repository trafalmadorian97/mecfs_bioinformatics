"""
Verify the plt.close fix: run the real SUSIE stackplot task repeatedly in one process
and confirm anonymous RSS returns to baseline after each plot, instead of climbing by
the retained figure's footprint (~1.8 GiB per plot) as it did before the fix.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.verify_stackplot_close \
    2>&1 | tee experiments/claude/genome_reference_harmonization/verify_stackplot_close.log
"""

import gc
import tempfile
from pathlib import Path

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr20_47_653_230_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES as G,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import make_wf

N_ITERS = 4


def _rss_anon_kib() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("RssAnon:"):
            return int(line.split()[1])
    raise RuntimeError("RssAnon missing")


def _fetch_from_store(task: Task):
    by_id = {dep.asset_id: dep for dep in task.deps}

    def fetch(asset_id: AssetId) -> Asset:
        meta = by_id[asset_id].meta
        path = DEFAULT_RUNNER.meta_to_path(meta)
        return DirectoryAsset(path) if isinstance(meta, DirMeta) else FileAsset(path)

    return fetch


def _gib(kib: float) -> float:
    return kib / 2**20


def main() -> None:
    task = G.susie_stackplot_task
    DEFAULT_RUNNER.run(list(task.deps))  # materialize inputs
    fetch = _fetch_from_store(task)
    gc.collect()
    print(f"baseline RSS {_gib(_rss_anon_kib()):.2f} GiB", flush=True)
    for i in range(N_ITERS):
        with tempfile.TemporaryDirectory() as scratch:
            task.execute(scratch_dir=Path(scratch), fetch=fetch, wf=make_wf())
        gc.collect()
        print(f"after plot {i + 1}: RSS {_gib(_rss_anon_kib()):.2f} GiB", flush=True)
    print("=== verify done ===", flush=True)


if __name__ == "__main__":
    main()
