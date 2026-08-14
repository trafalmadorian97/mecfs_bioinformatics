import hashlib
from collections.abc import Callable

from mecfs_bio.build_system.wf.object_store.base_object_store import (
    ObjectHead,
    ObjectStore,
)


class FakeObjectStore(ObjectStore):
    """In-memory ObjectStore for tests.

    head looks up the in-memory objects dict. upload_from_url does not fetch
    source_url; it records the call in .uploaded (in call order, duplicates
    included), records any on_progress callback in .on_progress_for keyed by uri,
    and registers a synthetic ObjectHead for uri so a subsequent head() finds it,
    mirroring the dedup contract a real store provides.
    """

    def __init__(self, objects: dict[str, ObjectHead] | None = None) -> None:
        self._objects: dict[str, ObjectHead] = dict(objects or {})
        self.uploaded: list[tuple[str, str]] = []
        self.on_progress_for: dict[str, Callable[[int], None] | None] = {}

    def head(self, uri: str) -> ObjectHead | None:
        return self._objects.get(uri)

    def upload_from_url(
        self,
        source_url: str,
        uri: str,
        on_progress: Callable[[int], None] | None = None,
    ) -> str:
        self.uploaded.append((source_url, uri))
        self.on_progress_for[uri] = on_progress
        sha256 = hashlib.sha256(f"{source_url}:{uri}".encode()).hexdigest()
        self._objects[uri] = ObjectHead(size_bytes=0, sha256=sha256)
        return sha256
