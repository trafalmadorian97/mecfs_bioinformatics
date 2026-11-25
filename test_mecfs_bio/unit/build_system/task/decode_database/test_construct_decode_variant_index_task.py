from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.decode_database.construct_decode_variant_index_task import (
    ConstructDecodeVariantIndexTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_INDEX_NAME_COL,
    DECODE_MISSING_ALLELE,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
)
from test_mecfs_bio.unit.build_system.task.decode_database.decode_fixtures import (
    FX_CLEAN_COL,
    FX_SPECIAL_COL,
    FX_TRUE_EAF_COL,
    FakeParallelRangeDownloader,
    FakeS3Access,
    membership_rows,
    synthetic_membership,
    synthetic_protein,
    synthetic_rows,
    write_clean_variant_info_asset,
    write_synthetic_decode_file,
)

_KEY = "final_somascan_smp/template.txt.gz"


def _meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=AssetId(asset_id),
        trait="t",
        project="p",
        sub_dir=Path("x"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def _build(tmp_path: Path, rows: pl.DataFrame, selected: pl.DataFrame) -> Path:
    source = tmp_path / "template.txt.gz"
    write_synthetic_decode_file(source, rows)
    membership_path = tmp_path / "membership.parquet"
    synthetic_membership(selected).write_parquet(membership_path)
    vi_path = tmp_path / "vi.parquet"
    write_clean_variant_info_asset(vi_path, rows)
    membership_task = FakeTask(meta=_meta("membership"))
    vi_task = FakeTask(meta=_meta("vi"))
    assets: dict[AssetId, Asset] = {
        membership_task.asset_id: FileAsset(membership_path),
        vi_task.asset_id: FileAsset(vi_path),
    }
    template = synthetic_protein(source, _KEY)
    task = ConstructDecodeVariantIndexTask.create(
        "index",
        template,
        membership_task=membership_task,
        clean_variant_info_task=vi_task,
    )
    wf = make_wf(
        s3_access=FakeS3Access(objects={template.ref: source}),
        parallel_range_downloader=FakeParallelRangeDownloader(),
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    asset = task.execute(scratch, fetch=lambda asset_id: assets[asset_id], wf=wf)
    assert isinstance(asset, FileAsset)
    return asset.path


@pytest.fixture(scope="module")
def built(tmp_path_factory: pytest.TempPathFactory) -> tuple[pl.DataFrame, Path]:
    rows = synthetic_rows(seed=7)
    path = _build(tmp_path_factory.mktemp("index"), rows, membership_rows(rows))
    return rows, path


def test_index_rows_are_the_clean_members_in_file_order(
    built: tuple[pl.DataFrame, Path],
) -> None:
    rows, path = built
    index = pl.read_parquet(path)
    members = membership_rows(rows)
    assert not members[FX_CLEAN_COL].all()
    assert index[DECODE_INDEX_NAME_COL].equals(
        members.filter(FX_CLEAN_COL)["Name"], check_names=False
    )
    bang = rows.filter(pl.col("otherAllele") == DECODE_MISSING_ALLELE)["Name"]
    assert not index[DECODE_INDEX_NAME_COL].is_in(bang.implode()).any()


def test_index_eaf_is_true_effect_allele_frequency(
    built: tuple[pl.DataFrame, Path],
) -> None:
    rows, path = built
    index = pl.read_parquet(path)
    truth = rows.select("Name", FX_TRUE_EAF_COL)
    joined = index.join(truth, left_on=DECODE_INDEX_NAME_COL, right_on="Name")
    assert index[GWASLAB_EFFECT_ALLELE_COL].null_count() == 0
    assert (joined[FX_TRUE_EAF_COL] > 0.5).any()
    assert (
        (joined[GWASLAB_EFFECT_ALLELE_FREQ_COL] - joined[FX_TRUE_EAF_COL]).abs() < 1e-6
    ).all()


def test_index_fails_when_eaf_undetermined(tmp_path: Path) -> None:
    rows = synthetic_rows(seed=7)
    ambiguous = rows.filter(pl.col(FX_SPECIAL_COL) & pl.col(FX_TRUE_EAF_COL).is_null())
    assert ambiguous.height == 1
    with pytest.raises(AssertionError):
        _build(tmp_path, rows, pl.concat([membership_rows(rows), ambiguous]))
