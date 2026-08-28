import gzip
from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genetic_map.parse_genetic_map_task import (
    GMAP_CM_COL,
    GMAP_POS_COL,
    GMAP_RATE_COL,
    ParseHg19GeneticMapTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL


def _write_map_gz(path: Path) -> None:
    # Eagle format: "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)", space-sep.
    # chr2 rows written before chr1 to prove the (CHR, POS) sort.
    lines = [
        "chr position COMBINED_rate(cM/Mb) Genetic_Map(cM)",
        "2 100 1.5 0.0",
        "1 721290 2.685807669 0.410292036939447",
        "1 55550 0 0",
    ]
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))


def test_parses_and_sorts_genetic_map(tmp_path: Path) -> None:
    raw = tmp_path / "genetic_map_hg19.txt.gz"
    _write_map_gz(raw)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    raw_task = FakeTask(
        ReferenceFileMeta(
            group="genetic_map",
            sub_group="hg19",
            sub_folder=PurePath("raw"),
            id=AssetId("genetic_map_raw"),
            extension=".txt.gz",
        )
    )
    task = ParseHg19GeneticMapTask.create(asset_id="genetic_map", raw_map_task=raw_task)

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "genetic_map_raw":
            return FileAsset(raw)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    df = pl.read_parquet(result.path)

    assert df.columns == [GWASLAB_CHROM_COL, GMAP_POS_COL, GMAP_RATE_COL, GMAP_CM_COL]
    assert df.schema[GWASLAB_CHROM_COL] == pl.Int64
    assert df.schema[GMAP_POS_COL] == pl.Int64
    assert df.schema[GMAP_RATE_COL] == pl.Float64
    # Sorted globally by (CHR, POS).
    assert df[GWASLAB_CHROM_COL].to_list() == [1, 1, 2]
    assert df[GMAP_POS_COL].to_list() == [55550, 721290, 100]
    # The recombination rate (cM/Mb) is carried through directly.
    assert abs(df[GMAP_RATE_COL][1] - 2.685807669) < 1e-9
