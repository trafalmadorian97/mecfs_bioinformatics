from pathlib import Path

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_S_LDSC,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from test_mecfs_bio.system.util import log_on_error


def test_run_s_lsdc(tmp_path: Path):
    """
    Test that we can run stratified linkage disequilibrium score regression on the DECODE ME data
    """
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"
    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run(DECODE_ME_S_LDSC.get_terminal_tasks())
        assert result is not None
