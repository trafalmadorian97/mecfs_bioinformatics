from abc import ABC, abstractmethod
from collections.abc import Callable

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
    def upload_from_url(
        self,
        source_url: str,
        uri: str,
        on_progress: Callable[[int], None] | None = None,
    ) -> str:
        """Stream source_url into uri and return the store's stored checksum.

        on_progress, if given, is invoked during the transfer with the number of
        bytes moved since the previous call (a delta, not a running total). It may be
        called from several threads concurrently, so an implementation must accept
        that and a callback must be thread-safe. The caller owns the callback and the
        meaning of progress (e.g. a percentage against a size it knows); the store
        only forwards byte deltas.
        """
        ...
