"""
End-to-end coverage of asset-store path remapping through SimpleRunner, which is the
level at which the feature is actually used.
"""

from pathlib import Path, PurePath
from typing import Any, Mapping, Sequence

import pytest
from structlog.testing import capture_logs

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    PathRemapRule,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.check_roots_available import (
    RemapRootUnavailableError,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.counting_task import CountingTask
from mecfs_bio.build_system.task.external_file_copy_task import ExternalFileCopyTask

# SimpleFileMeta assets live under this subtree of the store.
OTHER_FILES_PREFIX = PurePath("other_files")


def _warnings(logs: Sequence[Mapping[str, Any]]) -> list[Mapping[str, Any]]:
    """
    The warning-level events captured from structlog.  structlog is not routed through
    stdlib logging here, so pytest's caplog fixture sees nothing.
    """
    return [entry for entry in logs if entry.get("log_level") == "warning"]


def _remap_root(tmp_path: Path) -> Path:
    """
    Create the alternate root.  It has to exist before a run: an absent root is how a
    detached drive presents itself, and the runner refuses to start in that state.
    """
    remap_root = tmp_path / "remote_asset_store"
    remap_root.mkdir()
    return remap_root


def _external_file(tmp_path: Path) -> Path:
    external_dir = tmp_path / "external"
    external_dir.mkdir(parents=True, exist_ok=True)
    external_file = external_dir / "external_file.txt"
    external_file.write_text("abc123")
    return external_file


def test_remapped_asset_is_written_to_and_retrieved_from_the_remap_root(
    tmp_path: Path,
) -> None:
    """
    The asset must materialize under the remap root, and a second run must FIND it there
    rather than rebuild it.  Retrieval is the half that breaks if the mapper and
    get_asset_if_exists ever disagree about where an asset lives.
    """
    asset_root = tmp_path / "asset_store"
    remap_root = _remap_root(tmp_path)
    task = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("remapped_file")),
            external_path=_external_file(tmp_path),
        )
    )
    runner = SimpleRunner(
        info_store=tmp_path / "info_store.yaml",
        asset_root=asset_root,
        tracer=ImoHasher.with_xxhash_128(),
        path_remap=(PathRemapRule(root=remap_root, prefixes=(OTHER_FILES_PREFIX,)),),
    )

    store = runner.run([task])

    asset = store[task.asset_id]
    assert isinstance(asset, FileAsset)
    asset_path = asset.path
    assert asset_path == remap_root / OTHER_FILES_PREFIX / "remapped_file"
    assert asset_path.read_text() == "abc123"
    assert not (asset_root / OTHER_FILES_PREFIX).exists()
    assert task.run_count == 1

    runner.run([task])
    assert task.run_count == 1


def test_run_aborts_when_a_remap_root_is_missing(tmp_path: Path) -> None:
    """
    The detached-drive case.  It has to fail before anything is scheduled, because the
    alternative is not a crash but a quiet full rebuild onto the default root: assets on a
    drive that is not there are indistinguishable from assets that were never built.
    """
    task = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("remapped_file")),
            external_path=_external_file(tmp_path),
        )
    )
    runner = SimpleRunner(
        info_store=tmp_path / "info_store.yaml",
        asset_root=tmp_path / "asset_store",
        tracer=ImoHasher.with_xxhash_128(),
        path_remap=(
            PathRemapRule(
                root=tmp_path / "drive_not_attached",
                prefixes=(OTHER_FILES_PREFIX,),
            ),
        ),
    )

    with pytest.raises(RemapRootUnavailableError):
        runner.run([task])

    assert task.run_count == 0
    assert not (tmp_path / "asset_store").exists()


def test_run_warns_when_a_remapped_subtree_was_never_migrated(tmp_path: Path) -> None:
    """
    Changing the config without moving the data silently costs a full rebuild, so the
    runner has to say something.  The warning belongs on run() rather than at import of
    the default runner, so that merely importing the runner stays free of side effects.
    """
    asset_root = tmp_path / "asset_store"
    remap_root = _remap_root(tmp_path)
    stale_dir = asset_root / OTHER_FILES_PREFIX
    stale_dir.mkdir(parents=True)
    task = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("some_file")),
            external_path=_external_file(tmp_path),
        )
    )
    runner = SimpleRunner(
        info_store=tmp_path / "info_store.yaml",
        asset_root=asset_root,
        tracer=ImoHasher.with_xxhash_128(),
        path_remap=(PathRemapRule(root=remap_root, prefixes=(OTHER_FILES_PREFIX,)),),
    )

    with capture_logs() as logs:
        runner.run([task])

    assert len(_warnings(logs)) == 1


def test_run_is_silent_when_nothing_needs_migrating(tmp_path: Path) -> None:
    task = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("some_file")),
            external_path=_external_file(tmp_path),
        )
    )
    runner = SimpleRunner(
        info_store=tmp_path / "info_store.yaml",
        asset_root=tmp_path / "asset_store",
        tracer=ImoHasher.with_xxhash_128(),
        path_remap=(
            PathRemapRule(
                root=_remap_root(tmp_path),
                prefixes=(PurePath("reference_data"),),
            ),
        ),
    )

    with capture_logs() as logs:
        runner.run([task])

    assert _warnings(logs) == []


def test_unremapped_asset_stays_under_the_default_root(tmp_path: Path) -> None:
    asset_root = tmp_path / "asset_store"
    remap_root = _remap_root(tmp_path)
    task = CountingTask(
        ExternalFileCopyTask(
            meta=SimpleFileMeta(AssetId("local_file")),
            external_path=_external_file(tmp_path),
        )
    )
    runner = SimpleRunner(
        info_store=tmp_path / "info_store.yaml",
        asset_root=asset_root,
        tracer=ImoHasher.with_xxhash_128(),
        path_remap=(
            PathRemapRule(root=remap_root, prefixes=(PurePath("reference_data"),)),
        ),
    )

    store = runner.run([task])

    asset = store[task.asset_id]
    assert isinstance(asset, FileAsset)
    assert asset.path == asset_root / OTHER_FILES_PREFIX / "local_file"
    assert list(remap_root.iterdir()) == []
