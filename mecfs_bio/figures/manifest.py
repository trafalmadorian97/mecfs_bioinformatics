"""
Manifest for the figure system.

The manifest is a JSON file committed to git that maps figure paths (relative to
the figure directory) to the SHA-256 hash of the figure's content. It is the
source of truth for which figures are "live" in the project.

Figure blobs are stored content-addressed on the GitHub release: the asset name
on the release is the SHA-256 hash. Because storage is content-addressed and
append-only, two collaborators updating the same figure never overwrite each
other's blobs --- conflicts surface as ordinary git merge conflicts on the
manifest JSON, where one hash must be chosen.
"""

import hashlib
import json
from pathlib import Path

import structlog
from attrs import frozen

logger = structlog.get_logger()

MANIFEST_VERSION = 1
_HASH_FIELD = "figures"
_VERSION_FIELD = "version"


@frozen
class FigureManifest:
    """
    Mapping from figure path (relative to the figure directory) to the
    SHA-256 hex digest of file contents.

    Keys are ``Path`` objects in memory; they are serialised as POSIX
    strings in the on-disk JSON so the manifest is portable across
    operating systems.
    """

    figures: dict[Path, str]

    @classmethod
    def empty(cls) -> "FigureManifest":
        return cls(figures={})

    @classmethod
    def load(cls, manifest_path: Path) -> "FigureManifest":
        if not manifest_path.exists():
            return cls.empty()
        with open(manifest_path) as f:
            raw = json.load(f)
        version = raw.get(_VERSION_FIELD)
        if version != MANIFEST_VERSION:
            raise ValueError(
                f"Manifest at {manifest_path} has unsupported version {version}; "
                f"expected {MANIFEST_VERSION}"
            )
        figures = raw.get(_HASH_FIELD, {})
        assert isinstance(figures, dict)
        return cls(figures={Path(k): v for k, v in figures.items()})

    def save(self, manifest_path: Path) -> None:
        manifest_path.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            _VERSION_FIELD: MANIFEST_VERSION,
            _HASH_FIELD: {
                path.as_posix(): sha for path, sha in sorted(self.figures.items())
            },
        }
        with open(manifest_path, "w") as f:
            json.dump(payload, f, indent=2, sort_keys=False)
            f.write("\n")

    def with_entry(self, rel_path: Path, sha256: str) -> "FigureManifest":
        new = dict(self.figures)
        new[rel_path] = sha256
        return FigureManifest(figures=new)

    def without_entry(self, rel_path: Path) -> "FigureManifest":
        new = dict(self.figures)
        new.pop(rel_path, None)
        return FigureManifest(figures=new)

    def hashes(self) -> set[str]:
        return set(self.figures.values())


def sha256_of_file(path: Path, chunk_size: int = 1 << 16) -> str:
    hasher = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            chunk = f.read(chunk_size)
            if not chunk:
                break
            hasher.update(chunk)
    return hasher.hexdigest()


def scan_figure_dir(fig_dir: Path) -> FigureManifest:
    """
    Build a manifest by hashing every file under ``fig_dir`` recursively.
    Paths are stored as ``Path`` objects relative to ``fig_dir``.
    """
    figures: dict[Path, str] = {}
    for path in sorted(fig_dir.rglob("*")):
        if not path.is_file():
            continue
        rel = path.relative_to(fig_dir)
        figures[rel] = sha256_of_file(path)
    return FigureManifest(figures=figures)
