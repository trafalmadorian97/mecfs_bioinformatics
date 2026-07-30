import gzip
import tarfile
from pathlib import Path

import pytest

from mecfs_bio.build_system.task.ppp_database.protein_sample_size_task import (
    _first_regular_member_data_offset,
    extract_regenie_n_from_tar_head,
)

_REGENIE_HEADER = (
    "CHROM GENPOS ID ALLELE0 ALLELE1 A1FREQ INFO N TEST BETA SE CHISQ LOG10P EXTRA"
)


def _make_protein_tar(tmp_path: Path, n: int) -> Path:
    """A protein tar shaped like the real UKB-PPP upload: a directory member followed by
    per-chromosome gzipped regenie files, N constant across rows."""
    protein_dir = tmp_path / "PROTEIN"
    protein_dir.mkdir()
    for chrom in (1, 2):
        with gzip.open(protein_dir / f"discovery_chr{chrom}_x.gz", "wt") as handle:
            handle.write(_REGENIE_HEADER + "\n")
            for pos in range(1, 2001):  # enough rows that the head spans many lines
                handle.write(
                    f"{chrom} {pos} v:{pos} G A 0.2 0.9 {n} ADD 0.1 0.2 1 1 NA\n"
                )
    tar_path = tmp_path / "protein.tar"
    # ustar (not Python's default pax) mirrors the real UKB-PPP tars: a 512-byte
    # directory stub, then the first regular file's data at byte 1024.
    with tarfile.open(tar_path, "w", format=tarfile.USTAR_FORMAT) as tar:
        tar.add(protein_dir, arcname="PROTEIN")
    return tar_path


def test_first_regular_member_data_offset_skips_directory_stub(tmp_path: Path):
    head = _make_protein_tar(tmp_path, n=33709).read_bytes()[:65536]
    # Directory stub header at 0, data at 512; first regular file's data at 1024.
    assert _first_regular_member_data_offset(head) == 1024


def test_extract_regenie_n_from_tar_head(tmp_path: Path):
    head = _make_protein_tar(tmp_path, n=33709).read_bytes()[:65536]
    assert extract_regenie_n_from_tar_head(head) == 33709


def test_extract_regenie_n_from_tar_head_other_value(tmp_path: Path):
    head = _make_protein_tar(tmp_path, n=41069).read_bytes()[:65536]
    assert extract_regenie_n_from_tar_head(head) == 41069


def test_first_regular_member_data_offset_raises_without_member():
    with pytest.raises(ValueError, match="no regular-file member"):
        _first_regular_member_data_offset(b"\x00" * 512)
