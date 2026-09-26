import os
from pathlib import Path

import xxhash
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.simple_hasher import (
    HashConstructor,
    HashObject,
)


# Reads exactly sample_size bytes from file object f, erroring if EOF is reached before that.
def _read_sample_size(f, sample_size: int) -> bytes:
    data = b""
    while len(data) < sample_size:
        chunk = f.read(sample_size - len(data))
        if not chunk:
            raise EOFError("Could not read enough data")
        if len(chunk) == sample_size:
            return chunk
        data += chunk
    return data


@frozen(slots=True)
class ImoHasher(Tracer):
    """
    Based on : https://github.com/kalafut/py-imohash/blob/master/imohash/imohash.py
    Quick approximate hash that subsamples large files
    """

    hash_constructor: HashConstructor
    hash_name: str
    sample_threshold: int = 128 * 1024
    sample_size: int = 16 * 1024

    @classmethod
    def with_xxhash_32(cls):
        return cls(hash_constructor=xxhash.xxh32, hash_name="xxh32")

    @classmethod
    def with_xxhash_128(cls):
        return cls(hash_constructor=xxhash.xxh3_128, hash_name="xxh3_128")

    def __call__(self, a: Asset) -> str:
        if isinstance(a, FileAsset):
            file_hash = self.hash_constructor()
            hash, size = _process_file(
                a.path,
                file_hash,
                sample_size=self.sample_size,
                sample_threshold=self.sample_threshold,
            )
            return _finish(hash, size, hash_name=self.hash_name)

        if isinstance(a, DirectoryAsset):
            dir_hash = self.hash_constructor()
            total_size = 0
            for file_path in sorted(a.path.rglob("*")):
                if file_path.is_file():
                    dir_hash, inc_size = _process_file(
                        pth=file_path,
                        hash_object=dir_hash,
                        sample_size=self.sample_size,
                        sample_threshold=self.sample_threshold,
                    )
                    total_size += inc_size

            return _finish(dir_hash, total_size, hash_name=self.hash_name)

        raise ValueError(f"Unknown asset {a} of type{type(a)}")


def _process_file(
    pth: Path, hash_object: HashObject, sample_size: int, sample_threshold: int
) -> tuple[HashObject, int]:
    with open(pth, "rb") as f:
        f.seek(0, os.SEEK_END)
        size = f.tell()
        f.seek(0, os.SEEK_SET)
        if size < sample_threshold or sample_size < 1 or size < (4 * sample_size):
            data = f.read()
        else:
            data = _read_sample_size(f, sample_size)
            f.seek(size // 2)
            data += _read_sample_size(f, sample_size)
            f.seek(-sample_size, os.SEEK_END)
            data += _read_sample_size(f, sample_size)
    hash_object.update(data)
    return hash_object, size


def _finish(hash_object: HashObject, size: int, hash_name: str) -> str:
    return f"imo_{hash_name}:{hash_object.hexdigest()}_sz:{size}"
