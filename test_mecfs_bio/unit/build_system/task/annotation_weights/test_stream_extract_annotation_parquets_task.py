import io
import tarfile
from pathlib import Path, PurePath

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.annotation_weights.stream_extract_annotation_parquets_task import (
    LDSCORE_PARQUET_MEMBER_RE,
    StreamExtractAnnotationParquetsTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.polyfun_constants import (
    POLYFUN_A1_COL,
    POLYFUN_A2_COL,
    POLYFUN_BP_COL,
    POLYFUN_CHR_COL,
    POLYFUN_SNP_COL,
)


def _unused_fetch(asset_id: AssetId) -> Asset:
    raise AssertionError("this task has no dependencies; fetch must not be called")


def _member_bytes(frame: pl.DataFrame) -> bytes:
    buf = io.BytesIO()
    frame.write_parquet(buf)
    return buf.getvalue()


def _build_tarball(path: Path, members: dict[str, bytes]) -> None:
    with tarfile.open(path, mode="w:gz") as tar:
        for name, data in members.items():
            info = tarfile.TarInfo(name=name)
            info.size = len(data)
            tar.addfile(info, io.BytesIO(data))


def _meta() -> ReferenceDataDirectoryMeta:
    return ReferenceDataDirectoryMeta(
        group="polyfun",
        sub_group="annotations",
        sub_folder=PurePath("raw"),
        id=AssetId("annot_parquets"),
    )


def _annot_frame(chrom: int) -> pl.DataFrame:
    return pl.DataFrame(
        {
            POLYFUN_CHR_COL: [chrom, chrom],
            POLYFUN_BP_COL: [1, 2],
            POLYFUN_SNP_COL: ["rsA", "rsB"],
            POLYFUN_A1_COL: ["A", "G"],
            POLYFUN_A2_COL: ["C", "T"],
        }
    )


def test_extracts_only_annot_parquet_members(tmp_path: Path):
    tarball = tmp_path / "bundle.tar.gz"
    _build_tarball(
        tarball,
        {
            # arbitrary member-path prefix, interleaved with a decoy ld-score member
            "UKBB_LD/baselineLF2.2.UKB.1.annot.parquet": _member_bytes(_annot_frame(1)),
            "UKBB_LD/baselineLF2.2.UKB.1.l2.ldscore.parquet": _member_bytes(
                pl.DataFrame({"junk": [0]})
            ),
            "UKBB_LD/baselineLF2.2.UKB.2.annot.parquet": _member_bytes(_annot_frame(2)),
            "UKBB_LD/baselineLF2.2.UKB.1.l2.M": b"1\t2\t3\n",
        },
    )

    task = StreamExtractAnnotationParquetsTask(
        meta=_meta(),
        url="http://example.invalid/bundle.tar.gz",
        stream_opener=lambda _url: open(tarball, "rb"),
        required_chromosomes=frozenset({1, 2}),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=_unused_fetch, wf=make_wf())

    assert isinstance(result, DirectoryAsset)
    names = sorted(p.name for p in result.path.iterdir())
    assert names == [
        "baselineLF2.2.UKB.1.annot.parquet",
        "baselineLF2.2.UKB.2.annot.parquet",
    ]
    # the extracted member is a faithful copy of the original parquet
    got = pl.read_parquet(result.path / "baselineLF2.2.UKB.1.annot.parquet")
    assert got.columns == [
        POLYFUN_CHR_COL,
        POLYFUN_BP_COL,
        POLYFUN_SNP_COL,
        POLYFUN_A1_COL,
        POLYFUN_A2_COL,
    ]
    assert got[POLYFUN_A2_COL].to_list() == ["C", "T"]


def test_raises_when_a_required_chromosome_is_missing(tmp_path: Path):
    tarball = tmp_path / "bundle.tar.gz"
    _build_tarball(
        tarball,
        {"baselineLF2.2.UKB.1.annot.parquet": _member_bytes(_annot_frame(1))},
    )
    task = StreamExtractAnnotationParquetsTask(
        meta=_meta(),
        url="http://example.invalid/bundle.tar.gz",
        stream_opener=lambda _url: open(tarball, "rb"),
        required_chromosomes=frozenset({1, 2}),  # 2 is absent
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    with pytest.raises(ValueError):
        task.execute(scratch_dir=scratch, fetch=_unused_fetch, wf=make_wf())


def _ldscore_frame(chrom: int) -> pl.DataFrame:
    return pl.DataFrame(
        {
            POLYFUN_CHR_COL: [chrom, chrom],
            POLYFUN_BP_COL: [1, 2],
            POLYFUN_SNP_COL: ["rsA", "rsB"],
            POLYFUN_A1_COL: ["A", "G"],
            POLYFUN_A2_COL: ["C", "T"],
            "Coding_UCSC_common": [0.1, 0.2],
        }
    )


def test_extracts_only_ldscore_members(tmp_path: Path):
    tarball = tmp_path / "bundle.tar.gz"
    _build_tarball(
        tarball,
        {
            "UKBB_LD/baselineLF2.2.UKB.1.annot.parquet": _member_bytes(_annot_frame(1)),
            "UKBB_LD/baselineLF2.2.UKB.1.l2.ldscore.parquet": _member_bytes(
                _ldscore_frame(1)
            ),
            "UKBB_LD/baselineLF2.2.UKB.2.l2.ldscore.parquet": _member_bytes(
                _ldscore_frame(2)
            ),
            "UKBB_LD/baselineLF2.2.UKB.1.l2.M": b"1\t2\t3\n",
        },
    )
    task = StreamExtractAnnotationParquetsTask(
        meta=_meta(),
        url="http://example.invalid/bundle.tar.gz",
        stream_opener=lambda _url: open(tarball, "rb"),
        required_chromosomes=frozenset({1, 2}),
        member_pattern=LDSCORE_PARQUET_MEMBER_RE,
        dest_stem="baselineLF2.2.UKB.{chrom}.l2.ldscore",
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=_unused_fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    names = sorted(p.name for p in result.path.iterdir())
    assert names == [
        "baselineLF2.2.UKB.1.l2.ldscore.parquet",
        "baselineLF2.2.UKB.2.l2.ldscore.parquet",
    ]
