from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.annotation_weights.annotation_weights_comparison_task import (
    GAMMA_EXT_COL,
    GAMMA_INT_COL,
    AnnotationWeightsComparisonTask,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_RAW_COL,
    GAMMA_STANDARDIZED_COL,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf

_ANNOTATIONS = ["a", "b", "c"]
_FAMILIES = ["coding", "conserved", "coding"]


def _weights_table_task(
    tmp_path: Path, name: str, gamma_raw: list[float]
) -> tuple[FakeTask, FileAsset]:
    path = tmp_path / f"{name}.parquet"
    pl.DataFrame(
        {
            ANNOTATION_COL: _ANNOTATIONS,
            GAMMA_RAW_COL: gamma_raw,
            GAMMA_STANDARDIZED_COL: [0.0] * len(_ANNOTATIONS),
            FAMILY_COL: _FAMILIES,
        }
    ).write_parquet(path)
    task = FakeTask(
        ResultTableMeta(id=AssetId(name), trait="t", project="p", extension=".parquet")
    )
    return task, FileAsset(path)


def test_gammas_averaged_across_tables_and_scaled_to_max_magnitude(tmp_path: Path):
    ext, ext_asset = _weights_table_task(tmp_path, "ext", [1.0, -4.0, 2.0])
    odd, odd_asset = _weights_table_task(tmp_path, "odd", [1.0, 0.0, 2.0])
    even, even_asset = _weights_table_task(tmp_path, "even", [3.0, 1.0, 0.0])
    fetch_map: dict[str, Asset] = {
        "ext": ext_asset,
        "odd": odd_asset,
        "even": even_asset,
    }

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    task = AnnotationWeightsComparisonTask.create(
        asset_id="comparison",
        external_weights_tasks=(ext,),
        internal_weights_tasks=(odd, even),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    table = pl.read_parquet(result.path)

    by_annotation = {row[ANNOTATION_COL]: row for row in table.iter_rows(named=True)}
    # External (1, -4, 2) over max |gamma| 4.
    assert by_annotation["a"][GAMMA_EXT_COL] == pytest.approx(0.25)
    assert by_annotation["b"][GAMMA_EXT_COL] == pytest.approx(-1.0)
    assert by_annotation["c"][GAMMA_EXT_COL] == pytest.approx(0.5)
    # Internal: mean of odd and even is (2, 0.5, 1), over max |gamma| 2.
    assert by_annotation["a"][GAMMA_INT_COL] == pytest.approx(1.0)
    assert by_annotation["b"][GAMMA_INT_COL] == pytest.approx(0.25)
    assert by_annotation["c"][GAMMA_INT_COL] == pytest.approx(0.5)
    assert by_annotation["b"][FAMILY_COL] == "conserved"
    # c has the smallest larger-magnitude coefficient, so it sorts last.
    assert table[ANNOTATION_COL][-1] == "c"
