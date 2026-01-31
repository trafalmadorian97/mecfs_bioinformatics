import tempfile
from pathlib import Path

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.magma.decode_me_hba_magma_analysis import (
    DECODE_ME_HBA_MAGMA_TASKS,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner


def test_run_hba_magma():
    """
    Test that we can run MAGMA via the human brain atlas on the DECODE ME data
    """
    with tempfile.TemporaryDirectory() as tempdirname:
        tempdir = Path(tempdirname)
        info_store = tempdir / "info_store.yaml"
        asset_root = tempdir / "asset_store"
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(
            DECODE_ME_HBA_MAGMA_TASKS.terminal_tasks(),
        )
        assert result is not None
