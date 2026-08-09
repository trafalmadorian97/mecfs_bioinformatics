from abc import ABC, abstractmethod

from attrs import frozen


@frozen
class ObjectHead:
    """Metadata about an object already staged in a store, without its body."""

    size_bytes: int
    sha256: str | None


class ObjectStore(ABC):
    """Stages large files into a remote object store (e.g. S3) and checks presence.

    head lets a caller determine whether an object already exists at a uri
    (and, if so, its size and checksum) without downloading it, so callers can
    dedup a staging step. upload_from_url streams a remote source directly into
    the store, so the caller never has to hold the object body in memory or on
    local disk.
    """

    @abstractmethod
    def head(self, uri: str) -> ObjectHead | None: ...

    @abstractmethod
    def upload_from_url(self, source_url: str, uri: str) -> str: ...
