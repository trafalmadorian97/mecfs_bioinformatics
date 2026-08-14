from pathlib import Path

import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.sbayesrc.gctb_gwfm_constants import (
    GWFM_REFERENCE_BUNDLE,
    GWFM_REFERENCE_VERSION,
    MARKER_VERSION_KEY,
    gwfm_reference_prefix,
)
from mecfs_bio.build_system.task.sbayesrc.stage_gwfm_reference_task import (
    StageGwfmReferenceTask,
    make_upload_progress_logger,
)
from mecfs_bio.build_system.wf.base_wf import WF, make_wf
from mecfs_bio.build_system.wf.object_store.base_object_store import ObjectHead
from mecfs_bio.build_system.wf.object_store.fake_object_store import FakeObjectStore


def _uri(bucket: str, filename: str) -> str:
    return f"s3://{bucket}/{gwfm_reference_prefix(GWFM_REFERENCE_VERSION)}{filename}"


def _noop_fetch(asset_id: AssetId) -> Asset:
    raise AssertionError(
        "StageGwfmReferenceTask has no deps; fetch should not be called"
    )


@pytest.fixture
def wf_all_present() -> WF:
    objects = {
        _uri("mybucket", f.filename): ObjectHead(size_bytes=f.size_bytes, sha256=None)
        for f in GWFM_REFERENCE_BUNDLE
    }
    return make_wf(object_store=FakeObjectStore(objects=objects))


@pytest.fixture
def wf_one_missing() -> WF:
    objects = {
        _uri("mybucket", f.filename): ObjectHead(size_bytes=f.size_bytes, sha256=None)
        for f in GWFM_REFERENCE_BUNDLE
    }
    del objects[_uri("mybucket", GWFM_REFERENCE_BUNDLE[0].filename)]
    return make_wf(object_store=FakeObjectStore(objects=objects))


def test_skips_upload_when_bucket_already_populated(
    tmp_path: Path, wf_all_present: WF
) -> None:
    task = StageGwfmReferenceTask.create(bucket="mybucket")
    asset = task.execute(tmp_path, fetch=_noop_fetch, wf=wf_all_present)
    object_store = wf_all_present.object_store
    assert isinstance(object_store, FakeObjectStore)
    assert object_store.uploaded == []
    assert MARKER_VERSION_KEY in Path(asset.path).read_text()


def test_uploads_only_missing_files(tmp_path: Path, wf_one_missing: WF) -> None:
    task = StageGwfmReferenceTask.create(bucket="mybucket")
    task.execute(tmp_path, fetch=_noop_fetch, wf=wf_one_missing)
    object_store = wf_one_missing.object_store
    assert isinstance(object_store, FakeObjectStore)
    assert len(object_store.uploaded) == 1


def test_upload_passes_a_progress_callback_for_the_uploaded_file(
    tmp_path: Path, wf_one_missing: WF
) -> None:
    # The uploaded (missing) file gets a progress callback so the multi-hour transfer
    # is observable; files that are skipped never trigger an upload and so get none.
    task = StageGwfmReferenceTask.create(bucket="mybucket")
    task.execute(tmp_path, fetch=_noop_fetch, wf=wf_one_missing)
    object_store = wf_one_missing.object_store
    assert isinstance(object_store, FakeObjectStore)
    missing_uri = _uri("mybucket", GWFM_REFERENCE_BUNDLE[0].filename)
    assert object_store.on_progress_for[missing_uri] is not None


def test_progress_logger_throttles_to_band_boundaries_and_completion() -> None:
    calls: list[dict] = []
    on_progress = make_upload_progress_logger(
        filename="f.zip",
        total_bytes=1000,
        log=lambda _event, **kwargs: calls.append(kwargs),
    )
    # 4% stays inside the first 5% band -> no log yet.
    on_progress(40)
    assert calls == []
    # crossing to 8% emits exactly one line at the band boundary.
    on_progress(40)
    assert len(calls) == 1
    assert calls[-1]["percent"] == 8
    # a single delta that jumps straight to the end logs once at completion, not once
    # per crossed band.
    on_progress(920)
    assert len(calls) == 2
    assert calls[-1]["percent"] == 100


def test_progress_logger_logs_once_when_a_single_delta_completes_the_file() -> None:
    calls: list[dict] = []
    on_progress = make_upload_progress_logger(
        filename="f.zip",
        total_bytes=1000,
        log=lambda _event, **kwargs: calls.append(kwargs),
    )
    on_progress(1000)
    assert len(calls) == 1
    assert calls[0]["percent"] == 100
