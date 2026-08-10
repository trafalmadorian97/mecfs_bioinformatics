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
