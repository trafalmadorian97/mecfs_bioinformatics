"""Stream a polyfun baseline-LF tarball and keep only one kind of per-chromosome
parquet member, selected by an injected member pattern.

The polyfun bundle (baselineLF_v2.2.UKB.polyfun.tar.gz, ~30GB) contains the ~0.7GB
of allele-bearing per-chromosome annotation parquets
(baselineLF2.2.UKB.<chr>.annot.parquet) and the ~29GB of per-chromosome annotation
LD-score parquets (baselineLF2.2.UKB.<chr>.l2.ldscore.parquet). gzip is not
seekable, so we read the tarball sequentially from the URL in streaming mode
(tarfile r|gz) and copy out only the members matching member_pattern as they pass,
never storing the whole tarball. The output is a directory holding one extracted
parquet per chromosome, named dest_stem (with the chromosome filled in) plus
".parquet". The defaults select the annotation members.
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
LDSCORE_PARQUET_MEMBER_RE = re.compile(
    r"baselineLF2\.2\.UKB\.(\d+)\.l2\.ldscore\.parquet$"
)

StreamOpener = Callable[[str], BinaryIO]


def _default_stream_opener(url: str) -> BinaryIO:
    return urllib.request.urlopen(url)


@frozen
class StreamExtractAnnotationParquetsTask(Task):
    meta: Meta
    url: str
    stream_opener: StreamOpener = field(default=_default_stream_opener, eq=False)
    required_chromosomes: frozenset[int] = frozenset(range(1, 23))
    # Group 1 of the pattern must capture the chromosome number.
    member_pattern: re.Pattern[str] = ANNOT_PARQUET_MEMBER_RE
    # {chrom} is filled with the matched chromosome; ".parquet" is appended.
    dest_stem: str = "baselineLF2.2.UKB.{chrom}.annot"

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
                    match = self.member_pattern.search(member.name)
                    if match is None:
                        continue
                    chrom = int(match.group(1))
                    dest = scratch_dir / (
                        self.dest_stem.format(chrom=chrom) + ".parquet"
                    )
                    source = tar.extractfile(member)
                    assert source is not None
                    with open(dest, "wb") as out:
                        shutil.copyfileobj(source, out)
                    found[chrom] = dest
                    logger.info(
                        "extracted parquet member",
                        chromosome=chrom,
                        member=member.name,
                        size_bytes=member.size,
                    )
        missing = sorted(self.required_chromosomes - set(found))
        if missing:
            raise ValueError(
                f"tarball at {self.url} is missing members matching "
                f"{self.member_pattern.pattern} for "
                f"chromosomes {missing}; found {sorted(found)}"
            )
        return DirectoryAsset(scratch_dir)
