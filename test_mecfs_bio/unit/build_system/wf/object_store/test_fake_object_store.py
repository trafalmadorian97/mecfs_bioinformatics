from mecfs_bio.build_system.wf.object_store.base_object_store import ObjectHead
from mecfs_bio.build_system.wf.object_store.fake_object_store import FakeObjectStore


def test_fake_head_and_upload_roundtrip() -> None:
    store = FakeObjectStore(
        objects={"s3://b/present": ObjectHead(size_bytes=10, sha256="a" * 64)}
    )
    assert store.head("s3://b/absent") is None
    present = store.head("s3://b/present")
    assert present is not None
    assert present.size_bytes == 10
    returned = store.upload_from_url("https://x/file", "s3://b/absent")
    assert ("https://x/file", "s3://b/absent") in store.uploaded
    assert store.head("s3://b/absent") is not None
    assert isinstance(returned, str)
