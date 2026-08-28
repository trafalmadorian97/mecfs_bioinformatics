from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    BuildBaselineLFAnnotationParquetTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _write_members(directory: Path) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    # chr1: a plain variant; a genuinely multiallelic site (rs2 at BP200 with two
    # DIFFERENT allele pairs and different annotations -> both must survive); and
    # an ordering-duplicate at BP300 (same unordered {C,G}, identical annotations
    # -> must collapse to one).
    pl.DataFrame(
        {
            "CHR": [1, 1, 1, 1, 1],
            "SNP": ["rs1", "rs2", "rs2b", "rs4", "rs4"],
            "BP": [100, 200, 200, 300, 300],
            "A1": ["A", "G", "G", "C", "G"],
            "A2": ["C", "T", "A", "G", "C"],
            "annotA": [1.0, 0.0, 5.0, 2.0, 2.0],
            "annotB": [0.5, 0.25, 0.9, 0.1, 0.1],
        }
    ).write_parquet(directory / "baselineLF2.2.UKB.1.annot.parquet")
    # chr2 written as a separate member to prove global (CHR,BP) ordering.
    pl.DataFrame(
        {
            "CHR": [2],
            "SNP": ["rs3"],
            "BP": [50],
            "A1": ["A"],
            "A2": ["T"],
            "annotA": [1.0],
            "annotB": [0.75],
        }
    ).write_parquet(directory / "baselineLF2.2.UKB.2.annot.parquet")


def test_builds_sorted_allele_bearing_annotation_parquet(tmp_path: Path) -> None:
    members = tmp_path / "members"
    _write_members(members)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    members_task = FakeTask(
        ReferenceDataDirectoryMeta(
            group="polyfun",
            sub_group="annotations",
            sub_folder=PurePath("raw"),
            id=AssetId("annot_members"),
        )
    )
    task = BuildBaselineLFAnnotationParquetTask.create(
        asset_id="annot_parquet", annot_members_task=members_task
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "annot_members":
            return DirectoryAsset(members)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    df = pl.read_parquet(result.path)

    # rs4 ordering-duplicate collapses (5 rows -> 4 on chr1) + 1 on chr2 = 5.
    assert df.height == 5
    # A1/A2 retained; no CM; annotations float32.
    assert set(df.columns) == {"CHR", "BP", "SNP", "A1", "A2", "annotA", "annotB"}
    assert df.schema["annotA"] == pl.Float32
    assert df.schema["annotB"] == pl.Float32
    # Sorted globally by (CHR, BP).
    assert df["CHR"].to_list() == [1, 1, 1, 1, 2]
    assert df["BP"].to_list() == [100, 200, 200, 300, 50]
    # Multiallelic BP200 keeps BOTH alleles with their distinct annotations.
    bp200 = df.filter((pl.col("CHR") == 1) & (pl.col("BP") == 200)).sort("annotA")
    assert bp200.height == 2
    assert bp200["annotA"].to_list() == [0.0, 5.0]
    # Unique on (CHR, BP, unordered allele key).
    ak = df.with_columns(
        (pl.min_horizontal("A1", "A2") + "_" + pl.max_horizontal("A1", "A2")).alias(
            "ak"
        )
    )
    assert ak.select("CHR", "BP", "ak").n_unique() == df.height
