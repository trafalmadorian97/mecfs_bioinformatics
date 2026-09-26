import hashlib
from pathlib import Path
from typing import Any, Callable

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.base_tracer import (
    Tracer,
)

HashObject = Any
HashConstructor = Callable[[], HashObject]


@frozen(slots=True)
class SimpleHasher(Tracer):
    hash_constructor: HashConstructor
    chunk_size: int = 8192

    def __call__(self, a: Asset) -> str:
        if isinstance(a, FileAsset):
            file_hash = self.hash_constructor()
            with open(Path(a.path), "rb") as f:
                while chunk := f.read(self.chunk_size):
                    file_hash.update(chunk)
            return file_hash.hexdigest()
        if isinstance(a, DirectoryAsset):
            dir_hash = self.hash_constructor()
            for file_path in sorted(a.path.rglob("*")):
                if file_path.is_file():
                    with open(file_path, "rb") as f:
                        while chunk := f.read(self.chunk_size):
                            dir_hash.update(chunk)
            return dir_hash.hexdigest()

        raise ValueError(f"Unknown asset {a} of type{type(a)}")

    @classmethod
    def md5_hasher(cls) -> "SimpleHasher":
        return cls(hashlib.md5)
