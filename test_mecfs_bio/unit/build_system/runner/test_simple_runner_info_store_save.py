"""
SimpleRunner skips rewriting an existing info store when a run changes nothing, and still
rewrites it when a run does change a trace.  A YAML comment appended to the store between
runs survives only if the store is not rewritten.
"""

from pathlib import Path

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask

SENTINEL = "# sentinel: present only if the store was not rewritten\n"


def test_info_store_rewritten_only_when_a_run_changes_it(tmp_path: Path) -> None:
    external_file = tmp_path / "external.txt"
    external_file.write_text("abc123")
    info_store = tmp_path / "info_store.yaml"
    task = ExternalFileCopyTask(
        meta=SimpleFileMeta(AssetId("copied_file")),
        external_path=external_file,
    )
    runner = SimpleRunner(
        info_store=info_store,
        asset_root=tmp_path / "asset_store",
        tracer=ImoHasher.with_xxhash_128(),
    )

    runner.run([task])
    with open(info_store, "a") as f:
        f.write(SENTINEL)

    runner.run([task])
    assert SENTINEL in info_store.read_text()

    external_file.write_text("changed")
    runner.run([task], must_rebuild_transitive=[task])
    assert SENTINEL not in info_store.read_text()
