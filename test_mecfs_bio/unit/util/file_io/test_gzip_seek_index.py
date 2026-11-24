import gzip
import random
from pathlib import Path

import indexed_gzip

from mecfs_bio.util.file_io.gzip_seek_index import (
    export_seek_index,
    load_seek_index,
    seek_points,
)

KIB = 1024


def _synthetic_text(n_bytes: int, seed: int) -> bytes:
    rng = random.Random(seed)
    lines: list[str] = []
    size = 0
    while size < n_bytes:
        line = "\t".join(
            [
                f"chr{rng.randint(1, 22)}",
                str(rng.randint(1, 250_000_000)),
                f"{rng.random():.6f}",
                f"{rng.gauss(0, 1):.5f}",
            ]
        )
        lines.append(line)
        size += len(line) + 1
    return ("\n".join(lines) + "\n").encode()


def test_rapidgzip_export_imports_into_indexed_gzip(tmp_path: Path) -> None:
    data = _synthetic_text(2 * 1024 * KIB, seed=1)
    gz_path = tmp_path / "data.txt.gz"
    gz_path.write_bytes(gzip.compress(data, compresslevel=1))
    index_path = tmp_path / "index.gzidx.zst"

    export_seek_index(gz_path, index_path, parallelization=4, chunk_size=32 * KIB)
    blob = index_path.read_bytes()

    with open(gz_path, "rb") as source:
        points = seek_points(blob, source)
    assert len(points) > 10

    rng = random.Random(2)
    with open(gz_path, "rb") as source:
        f = indexed_gzip.IndexedGzipFile(fileobj=source, auto_build=False)
        f.import_index(fileobj=load_seek_index(blob))
        for _ in range(50):
            offset = rng.randrange(len(data))
            length = rng.randint(1, 64 * KIB)
            f.seek(offset)
            assert f.read(length) == data[offset : offset + length]
        f.close()
