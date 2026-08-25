import gzip
import tarfile
from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    BuildBaselineLFAnnotationParquetTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _write_annot_gz(path: Path, rows: list[dict]) -> None:
    header = ["CHR", "BP", "SNP", "CM", "annotA", "annotB"]
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(r[c]) for c in header))
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))


def _build_fake_tarball(tmp_path: Path) -> Path:
    src = tmp_path / "src" / "baselineLF_v2.2.UKB"
    src.mkdir(parents=True)
    # chr2 out of order relative to chr1 to prove global (CHR,BP) sort;
    # a duplicate rsid on chr1 to prove dedup.
    _write_annot_gz(
        src / "baselineLF2.2.UKB.1.annot.gz",
        [
            {"CHR": 1, "BP": 100, "SNP": "rs1", "CM": 0.1, "annotA": 1, "annotB": 0.5},
            {"CHR": 1, "BP": 200, "SNP": "rs2", "CM": 0.2, "annotA": 0, "annotB": 0.25},
            {"CHR": 1, "BP": 200, "SNP": "rs2", "CM": 0.2, "annotA": 0, "annotB": 0.25},
        ],
    )
    _write_annot_gz(
        src / "baselineLF2.2.UKB.2.annot.gz",
        [
            {"CHR": 2, "BP": 50, "SNP": "rs3", "CM": 0.3, "annotA": 1, "annotB": 0.75},
        ],
    )
    # decoy non-annotation member that must be ignored
    (src / "baselineLF2.2.UKB.1.l2.ldscore.gz").write_bytes(gzip.compress(b"junk\n"))
    tarball = tmp_path / "baselineLF_v2.2.UKB.tar.gz"
    with tarfile.open(tarball, "w:gz") as tar:
        tar.add(src, arcname="baselineLF_v2.2.UKB")
    return tarball


def test_builds_sorted_deduped_annotation_parquet(tmp_path: Path) -> None:
    tarball = _build_fake_tarball(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    tarball_task = FakeTask(
        ReferenceFileMeta(
            group="polyfun",
            sub_group="annotations",
            sub_folder=PurePath("raw"),
            id=AssetId("annot_tarball"),
            extension=".tar.gz",
        )
    )
    task = BuildBaselineLFAnnotationParquetTask.create(
        asset_id="annot_parquet", tarball_task=tarball_task
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "annot_tarball":
            return FileAsset(tarball)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    df = pl.read_parquet(result.path)

    # dedup: rs2 appears once -> 3 unique SNPs total
    assert df.height == 3
    assert df["SNP"].to_list() == ["rs1", "rs2", "rs3"]  # (CHR,BP) sorted
    assert df["CHR"].to_list() == [1, 1, 2]
    # annotation columns present and float32
    assert df.schema["annotA"] == pl.Float32
    assert df.schema["annotB"] == pl.Float32
    assert set(df.columns) == {"CHR", "BP", "SNP", "CM", "annotA", "annotB"}
