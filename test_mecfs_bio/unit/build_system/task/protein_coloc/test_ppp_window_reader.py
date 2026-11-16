import gzip
import random
import tarfile
from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_reader import (
    PppWindowReader,
)
from mecfs_bio.build_system.wf.base_wf import WF, make_wf
from mecfs_bio.build_system.wf.synapse_downloader import SynapseDownloader
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.ppp_database_constants import Oid, SynID
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval
from mecfs_bio.util.file_io.ranged_tar import RangeReader

_REGENIE_HEADER = (
    "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA"
)

_PROTEIN = PppProteinFile(
    gene="GENE1",
    uniprot="P00001",
    oid=Oid("OID1"),
    version="v1",
    panel="Oncology",
    synid=SynID("syn1"),
)


@frozen(slots=True)
class _ByteRange:
    start: int
    end: int

    def overlap_length(self, other: "_ByteRange") -> int:
        return max(0, min(self.end, other.end) - max(self.start, other.start))


class _RecordingRangeReader:
    def __init__(self, data: bytes, reads: list[_ByteRange]) -> None:
        self.data = data
        self.reads = reads

    def read_range(self, start: int, length: int) -> bytes:
        self.reads.append(_ByteRange(start=start, end=start + length))
        return self.data[start : start + length]


class _TarServingDownloader(SynapseDownloader):
    """Serves byte ranges of a prepared tar, records them, and refuses full downloads."""

    def __init__(self, tar_path: Path) -> None:
        self.tar_path = tar_path
        self.reads: list[_ByteRange] = []

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        raise AssertionError("the fetcher must not download the whole tar")

    def open_range_reader(self, synid: str) -> RangeReader:
        return _RecordingRangeReader(self.tar_path.read_bytes(), self.reads)


def _write_chr_gz(path: Path, rows: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(_REGENIE_HEADER + "\n")
        handle.write("\n".join(rows) + "\n")


# Long enough that GNU tar stores member names in separate long-name entries, as in the
# real tars of complexes such as AMY1A_AMY1B_AMY1C.
_LONG_DIR_NAME = "GENE1_GENE2_GENE3_P00001_P00002_P00003_OID1_v1_Oncology_II"


def _filler_rows(chrom: int) -> list[str]:
    """Rows far from every test window, numerous and varied enough that the chromosome's
    compressed file is much larger than one tar header read."""
    rng = random.Random(0)
    return [
        f"{chrom} {10_000 + i} f{i} C T {rng.random():.6f} 0.9 1000 ADD "
        f"{rng.random():.6f} {rng.random():.6f} 1 1 NA"
        for i in range(20_000)
    ]


def _make_tar(tmp_path: Path) -> Path:
    member_dir = tmp_path / "build" / _LONG_DIR_NAME
    member_dir.mkdir(parents=True)
    _write_chr_gz(
        member_dir / "discovery_chr1_GENE1_P00001_OID1_v1_Oncology.gz",
        [
            "1 100 a G A 0.2 0.9 1000 ADD 0.5 0.1 1 1 NA",
            "1 200 b C T 0.3 0.9 1000 ADD 0.3 0.2 1 1 NA",
            "1 5000 c C T 0.3 0.9 1000 ADD 0.3 0.2 1 1 NA",
        ],
    )
    _write_chr_gz(
        member_dir / "discovery_chr2_GENE1_P00001_OID1_v1_Oncology.gz",
        ["2 300 d T A 0.1 0.9 1000 ADD 0.7 0.3 1 1 NA", *_filler_rows(chrom=2)],
    )
    _write_chr_gz(
        member_dir / "discovery_chrX_GENE1_P00001_OID1_v1_Oncology.gz",
        ["23 400 e T A 0.1 0.9 1000 ADD 0.7 0.3 1 1 NA"],
    )
    tar_path = tmp_path / "source.tar"
    with tarfile.open(tar_path, "w", format=tarfile.GNU_FORMAT) as tar:
        tar.add(member_dir, arcname=member_dir.name)
    return tar_path


def _setup(tmp_path: Path) -> tuple[_TarServingDownloader, WF, Path]:
    downloader = _TarServingDownloader(_make_tar(tmp_path))
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    return downloader, make_wf(synapse_downloader=downloader), scratch


_WINDOWS = [
    GenomicInterval(chrom=1, start=50, end=250),
    GenomicInterval(chrom=2, start=0, end=1000),
]


def test_fetch_returns_only_window_rows_in_common_schema(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    frames = PppWindowReader(proteins=(_PROTEIN,)).fetch_protein(
        "OID1", _WINDOWS, scratch, wf
    )
    assert set(frames) == set(_WINDOWS)
    l1 = frames[_WINDOWS[0]]
    assert l1.schema == WINDOW_SCHEMA
    assert l1[GWASLAB_POS_COL].to_list() == [100, 200]
    first = l1.row(0, named=True)
    assert first[GWASLAB_EFFECT_ALLELE_COL] == "A"
    assert first[GWASLAB_NON_EFFECT_ALLELE_COL] == "G"
    assert first[GWASLAB_EFFECT_ALLELE_FREQ_COL] == 0.2
    assert first[GWASLAB_BETA_COL] == 0.5
    assert frames[_WINDOWS[1]][GWASLAB_POS_COL].to_list() == [300]


def test_fetch_reads_chromosome_x_as_23(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    window = GenomicInterval(chrom=23, start=0, end=1000)
    frames = PppWindowReader(proteins=(_PROTEIN,)).fetch_protein(
        "OID1", [window], scratch, wf
    )
    assert frames[window][GWASLAB_POS_COL].to_list() == [400]


def _member_data_range(tar_path: Path, chrom_tag: str) -> _ByteRange:
    with tarfile.open(tar_path) as tar:
        (member,) = [m for m in tar.getmembers() if chrom_tag in m.name]
    return _ByteRange(start=member.offset_data, end=member.offset_data + member.size)


def test_fetch_does_not_read_an_unneeded_chromosome_file(tmp_path: Path) -> None:
    downloader, wf, scratch = _setup(tmp_path)
    PppWindowReader(proteins=(_PROTEIN,)).fetch_protein(
        "OID1", [GenomicInterval(chrom=1, start=50, end=250)], scratch, wf
    )
    unneeded = _member_data_range(downloader.tar_path, "_chr2_")
    bytes_read = sum(read.overlap_length(unneeded) for read in downloader.reads)
    # Reading the next tar header may touch the start of chr2's data, never all of it.
    assert bytes_read < (unneeded.end - unneeded.start) / 2


def test_window_with_no_rows_returns_empty_frame_with_schema(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    window = GenomicInterval(chrom=1, start=1000, end=2000)
    frames = PppWindowReader(proteins=(_PROTEIN,)).fetch_protein(
        "OID1", [window], scratch, wf
    )
    assert frames[window].height == 0
    assert frames[window].schema == WINDOW_SCHEMA
    assert isinstance(frames[window], pl.DataFrame)
