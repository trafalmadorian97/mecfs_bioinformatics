import gzip
import tarfile
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    PppProteinFile,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import (
    LocusWindow,
    PppWindowFetcher,
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


class _TarCopyingDownloader(SynapseDownloader):
    """Serves a prepared tar for every request and counts downloads."""

    def __init__(self, tar_path: Path) -> None:
        self.tar_path = tar_path
        self.n_calls = 0

    def download(self, synid: str, dest_dir: Path, expected_name: str) -> Path:
        self.n_calls += 1
        out = dest_dir / expected_name
        out.write_bytes(self.tar_path.read_bytes())
        return out

    def fetch_file_head(self, synid: str, n_bytes: int) -> bytes:
        return self.tar_path.read_bytes()[:n_bytes]


def _write_chr_gz(path: Path, rows: list[str]) -> None:
    with gzip.open(path, "wt") as handle:
        handle.write(_REGENIE_HEADER + "\n")
        handle.write("\n".join(rows) + "\n")


def _make_tar(tmp_path: Path) -> Path:
    member_dir = tmp_path / "build" / "GENE1_P00001_OID1_v1_Oncology"
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
        ["2 300 d T A 0.1 0.9 1000 ADD 0.7 0.3 1 1 NA"],
    )
    _write_chr_gz(
        member_dir / "discovery_chrX_GENE1_P00001_OID1_v1_Oncology.gz",
        ["23 400 e T A 0.1 0.9 1000 ADD 0.7 0.3 1 1 NA"],
    )
    tar_path = tmp_path / "source.tar"
    with tarfile.open(tar_path, "w") as tar:
        tar.add(member_dir, arcname=member_dir.name)
    return tar_path


def _setup(tmp_path: Path) -> tuple[_TarCopyingDownloader, WF, Path]:
    downloader = _TarCopyingDownloader(_make_tar(tmp_path))
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    return downloader, make_wf(synapse_downloader=downloader), scratch


_WINDOWS = [
    LocusWindow(locus_id="L1", chrom=1, start=50, end=250),
    LocusWindow(locus_id="L2", chrom=2, start=0, end=1000),
]


def test_fetch_returns_only_window_rows_in_common_schema(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    frames = PppWindowFetcher(proteins=(_PROTEIN,)).fetch("OID1", _WINDOWS, scratch, wf)
    assert set(frames) == {"L1", "L2"}
    l1 = frames["L1"]
    assert l1.schema == WINDOW_SCHEMA
    assert l1[GWASLAB_POS_COL].to_list() == [100, 200]
    first = l1.row(0, named=True)
    assert first[GWASLAB_EFFECT_ALLELE_COL] == "A"
    assert first[GWASLAB_NON_EFFECT_ALLELE_COL] == "G"
    assert first[GWASLAB_EFFECT_ALLELE_FREQ_COL] == 0.2
    assert first[GWASLAB_BETA_COL] == 0.5
    assert frames["L2"][GWASLAB_POS_COL].to_list() == [300]


def test_fetch_reads_chromosome_x_as_23(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    frames = PppWindowFetcher(proteins=(_PROTEIN,)).fetch(
        "OID1", [LocusWindow(locus_id="LX", chrom=23, start=0, end=1000)], scratch, wf
    )
    assert frames["LX"][GWASLAB_POS_COL].to_list() == [400]


def test_fetch_downloads_once_for_windows_on_two_chromosomes(tmp_path: Path) -> None:
    downloader, wf, scratch = _setup(tmp_path)
    PppWindowFetcher(proteins=(_PROTEIN,)).fetch("OID1", _WINDOWS, scratch, wf)
    assert downloader.n_calls == 1


def test_fetch_deletes_download_and_extraction(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    PppWindowFetcher(proteins=(_PROTEIN,)).fetch("OID1", _WINDOWS, scratch, wf)
    assert not any(scratch.rglob("*.tar"))
    assert not any(scratch.rglob("*.gz"))


def test_window_with_no_rows_returns_empty_frame_with_schema(tmp_path: Path) -> None:
    _, wf, scratch = _setup(tmp_path)
    frames = PppWindowFetcher(proteins=(_PROTEIN,)).fetch(
        "OID1", [LocusWindow(locus_id="E", chrom=1, start=1000, end=2000)], scratch, wf
    )
    assert frames["E"].height == 0
    assert frames["E"].schema == WINDOW_SCHEMA
    assert isinstance(frames["E"], pl.DataFrame)
