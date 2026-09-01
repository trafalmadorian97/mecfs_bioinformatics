"""Stream a polyfun baseline-LF tarball and keep only the per-chromosome
baselineLF2.2.UKB.<chr>.annot.parquet members.

The polyfun bundle (baselineLF_v2.2.UKB.polyfun.tar.gz, ~30GB) contains, besides
the ~0.7GB of allele-bearing annotation parquets we want, the ~29GB of
per-chromosome LD-score parquets that we do not need. gzip is not seekable, so
we read the tarball sequentially from the URL in streaming mode (tarfile r|gz)
and copy out only the matching members as they pass, never storing the whole
tarball. The output is a directory of the raw per-chromosome annotation
parquets, one per chromosome.
"""

import re
import shutil
import tarfile
import urllib.request
from collections.abc import Callable
from pathlib import Path
from typing import BinaryIO

import structlog
from attrs import field, frozen

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger(__name__)

ANNOT_PARQUET_MEMBER_RE = re.compile(r"baselineLF2\.2\.UKB\.(\d+)\.annot\.parquet$")

StreamOpener = Callable[[str], BinaryIO]


def _default_stream_opener(url: str) -> BinaryIO:
    return urllib.request.urlopen(url)


@frozen
class StreamExtractAnnotationParquetsTask(Task):
    meta: Meta
    url: str
    stream_opener: StreamOpener = field(default=_default_stream_opener, eq=False)
    required_chromosomes: frozenset[int] = frozenset(range(1, 23))

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> DirectoryAsset:
        found: dict[int, Path] = {}
        with self.stream_opener(self.url) as raw:
            with tarfile.open(fileobj=raw, mode="r|gz") as tar:
                for member in tar:
                    if not member.isfile():
                        continue
                    match = ANNOT_PARQUET_MEMBER_RE.search(member.name)
                    if match is None:
                        continue
                    chrom = int(match.group(1))
                    dest = scratch_dir / f"baselineLF2.2.UKB.{chrom}.annot.parquet"
                    source = tar.extractfile(member)
                    assert source is not None
                    with open(dest, "wb") as out:
                        shutil.copyfileobj(source, out)
                    found[chrom] = dest
                    logger.info(
                        "extracted annotation parquet",
                        chromosome=chrom,
                        member=member.name,
                        size_bytes=member.size,
                    )
        missing = sorted(self.required_chromosomes - set(found))
        if missing:
            raise ValueError(
                f"tarball at {self.url} is missing annot.parquet members for "
                f"chromosomes {missing}; found {sorted(found)}"
            )
        return DirectoryAsset(scratch_dir)
