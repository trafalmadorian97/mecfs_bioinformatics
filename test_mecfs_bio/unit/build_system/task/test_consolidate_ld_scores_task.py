"""
Simple LLM-implemented test
"""

import gzip
from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    ConsolidateLDScoresTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import SimpleWF


def _write_ld_score_gz(path: Path, rows: list[tuple[int, str, int, float]]):
    lines = ["CHR\tSNP\tBP\tL2\n"]
    for chrom, snp, bp, l2 in rows:
        lines.append(f"{chrom}\t{snp}\t{bp}\t{l2}\n")
    with gzip.open(path, "wt") as f:
        f.writelines(lines)


def test_consolidate_ld_scores(tmp_path: Path):
    source_dir = tmp_path / "ld_scores"
    source_dir.mkdir()

    _write_ld_score_gz(
        source_dir / "LDscore.1.l2.ldscore.gz",
        [
            (1, "rs111", 200, 80.5),
            (1, "rs112", 100, 90.3),
        ],
    )
    _write_ld_score_gz(
        source_dir / "LDscore.2.l2.ldscore.gz",
        [
            (2, "rs221", 300, 70.1),
            (2, "rs222", 150, 60.2),
        ],
    )

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

    result = tsk.execute(scratch_dir=scratch_dir, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)

    df = pl.read_parquet(result.path)
    assert len(df) == 4
    # Should be sorted by CHR then BP
    assert df["CHR"].to_list() == [1, 1, 2, 2]
    assert df["BP"].to_list() == [100, 200, 150, 300]
