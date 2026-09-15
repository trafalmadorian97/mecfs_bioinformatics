"""Task-level test of IndexedFastaTask, exercising the memory-mapped reference gather."""

import gzip
from pathlib import Path, PurePath

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.indexed_fasta_task import (
    IndexedFastaTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf

_LINE_WIDTH = 7
# 1-based: position 1 is A; 31-40 is a run of ten T.
_CHR1 = "ACGTTGGGGGCATGCAATTCGATCGATCGATTTTTTTTTTACACACACACGTCAGTCAGT"
_CHRX = "ggggaaaacccc"  # lowercase (soft-masked) bases must still match
_GZ_ID = "fasta_gz"


def _write_gzipped_fasta(path: Path) -> None:
    lines: list[str] = []
    for name, sequence in {"chr1": _CHR1, "chrX": _CHRX, "chrUn_gl000220": "ACGT"}.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[start : start + _LINE_WIDTH]
            for start in range(0, len(sequence), _LINE_WIDTH)
        )
    with gzip.open(path, "wt") as out:
        out.write("\n".join(lines) + "\n")


def _build_indexed_fasta(tmp_path: Path) -> IndexedFasta:
    gz_path = tmp_path / "genome.fa.gz"
    _write_gzipped_fasta(gz_path)
    source = FakeTask(
        ReferenceFileMeta(
            group="genome_sequence",
            sub_group="synthetic",
            sub_folder=PurePath("raw"),
            extension=".fa.gz",
            id=AssetId(_GZ_ID),
        )
    )
    task = IndexedFastaTask.create(fasta_gz_task=source, asset_id="indexed")

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == _GZ_ID
        return FileAsset(gz_path)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    return IndexedFasta.open(result.path)


def test_indexed_fasta_maps_main_contigs_to_gwaslab_codes(tmp_path: Path) -> None:
    fasta = _build_indexed_fasta(tmp_path)
    assert set(fasta.entries) == {1, 23}


@pytest.mark.parametrize("max_gather_bytes", [1, 1_000_000])
def test_reference_matches_across_line_wraps(tmp_path: Path, max_gather_bytes: int) -> None:
    fasta = _build_indexed_fasta(tmp_path)
    positions = np.array([1, 5, 5, 31, 31, 58, 60], dtype=np.int64)
    alleles = pl.Series(["A", "TG", "TC", "T" * 10, "T" * 11, "AGT", "TA"])
    matched = reference_matches(
        fasta,
        chrom=1,
        positions=positions,
        alleles=alleles,
        max_gather_bytes=max_gather_bytes,
    )
    # TA at 60 runs past the contig end, so it cannot match.
    assert matched.tolist() == [True, True, False, True, False, True, False]


def test_reference_matches_is_case_insensitive(tmp_path: Path) -> None:
    fasta = _build_indexed_fasta(tmp_path)
    matched = reference_matches(
        fasta,
        chrom=23,
        positions=np.array([1, 5], dtype=np.int64),
        alleles=pl.Series(["GGGG", "aaaac"]),
    )
    assert matched.tolist() == [True, True]
