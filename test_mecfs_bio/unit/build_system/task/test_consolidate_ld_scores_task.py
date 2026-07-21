"""
Simple LLM-implemented test
"""

import gzip
from pathlib import Path, PurePath

import narwhals
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_M_5_50_COL,
    ConsolidateLDScoresTask,
    total_m_5_50,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _write_ld_score_gz(path: Path, rows: list[tuple[int, str, int, float]]):
    lines = ["CHR\tSNP\tBP\tL2\n"]
    for chrom, snp, bp, l2 in rows:
        lines.append(f"{chrom}\t{snp}\t{bp}\t{l2}\n")
    with gzip.open(path, "wt") as f:
        f.writelines(lines)


def _write_ld_scores_with_counts(source_dir: Path):
    """One chromosome's LD scores plus the sibling M_5_50 file the task pairs them with."""
    _write_ld_score_gz(
        source_dir / "LDscore.1.l2.ldscore.gz",
        [
            (1, "rs111", 200, 80.5),
            (1, "rs112", 100, 90.3),
        ],
    )
    (source_dir / "LDscore.1.l2.M_5_50").write_text("1200\n")
    _write_ld_score_gz(
        source_dir / "LDscore.2.l2.ldscore.gz",
        [
            (2, "rs221", 300, 70.1),
            (2, "rs222", 150, 60.2),
        ],
    )
    (source_dir / "LDscore.2.l2.M_5_50").write_text("800\n")


def test_consolidate_ld_scores(tmp_path: Path):
    source_dir = tmp_path / "ld_scores"
    source_dir.mkdir()
    _write_ld_scores_with_counts(source_dir)

    source_meta = ReferenceDataDirectoryMeta(
        group="reference_data",
        sub_group="ld_scores",
        sub_folder=PurePath("extracted"),
        id=AssetId("dummy_extracted"),
    )
    fake_dep = FakeTask(meta=source_meta)

    tsk = ConsolidateLDScoresTask.create(
        asset_id="consolidated_ld",
        extracted_ld_score_task=fake_dep,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return DirectoryAsset(source_dir)

    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)

    df = pl.read_parquet(result.path)
    assert len(df) == 4
    # Should be sorted by CHR then BP
    assert df["CHR"].to_list() == [1, 1, 2, 2]
    assert df["BP"].to_list() == [100, 200, 150, 300]
    # Each chromosome's common-variant count is denormalized onto its own rows.
    assert df[LD_SCORE_M_5_50_COL].to_list() == [1200.0, 1200.0, 800.0, 800.0]
    assert total_m_5_50(narwhals.from_native(df, eager_only=True)) == 2000.0


def test_consolidate_ld_scores_requires_matching_variant_count(tmp_path: Path):
    """A chromosome whose M_5_50 file is missing must fail rather than silently contribute
    variants that the genome-wide count does not describe."""
    source_dir = tmp_path / "ld_scores"
    source_dir.mkdir()
    _write_ld_scores_with_counts(source_dir)
    (source_dir / "LDscore.2.l2.M_5_50").unlink()

    tsk = ConsolidateLDScoresTask.create(
        asset_id="consolidated_ld",
        extracted_ld_score_task=FakeTask(
            meta=ReferenceDataDirectoryMeta(
                group="reference_data",
                sub_group="ld_scores",
                sub_folder=PurePath("extracted"),
                id=AssetId("dummy_extracted"),
            )
        ),
    )
    scratch_dir = tmp_path / "scratch"
    scratch_dir.mkdir()
    with pytest.raises(AssertionError):
        tsk.execute(
            scratch_dir=scratch_dir,
            fetch=lambda asset_id: DirectoryAsset(source_dir),
            wf=make_wf(),
        )
