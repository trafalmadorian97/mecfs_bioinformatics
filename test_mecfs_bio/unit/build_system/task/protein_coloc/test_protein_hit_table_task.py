from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    BuildSlimProteinParquetTask,
    PppProteinFile,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    PppProteinHitTableTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ppp_database_constants import Oid, SynID
from mecfs_bio.constants.protein_coloc_constants import (
    HIT_CHR_COL,
    HIT_GENE_COL,
    HIT_LEAD_MLOG10P_COL,
    HIT_LEAD_POS_COL,
    HIT_N_SIG_VARIANTS_COL,
    HIT_PROTEIN_ID_COL,
    HIT_SEG_END_COL,
    HIT_SEG_START_COL,
)

_SE = 0.01
_NULL_Z = 0.5
# z = 6 gives p ~ 2e-9 < 5e-8; z = 5 gives p ~ 5.7e-7, significant only at 1e-5.
_HIT_Z = 6.0
_WEAK_Z = 5.0


def _protein(n: int) -> PppProteinFile:
    return PppProteinFile(
        gene=f"G{n}",
        uniprot=f"P{n}",
        oid=Oid(f"OID{n}"),
        version="v1",
        panel="Oncology",
        synid=SynID(f"syn{n}"),
    )


def _meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=AssetId(asset_id),
        trait="synthetic",
        project="test",
        sub_dir=PurePath("analysis"),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def _run(
    tmp_path: Path,
    variants: list[tuple[int, int]],
    z_by_protein: list[list[float]],
    p_threshold: float = 5e-8,
    merge_gap_bp: int = 250_000,
) -> pl.DataFrame:
    """variants: (CHR, POS) in index order; z_by_protein: one z per variant per protein."""
    tmp_path.mkdir(parents=True, exist_ok=True)
    assets: dict[AssetId, Asset] = {}
    index_task = FakeTask(meta=_meta("index"))
    index_path = tmp_path / "index.parquet"
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [c for c, _ in variants],
            GWASLAB_POS_COL: [p for _, p in variants],
        },
        schema={GWASLAB_CHROM_COL: pl.Int32, GWASLAB_POS_COL: pl.Int32},
    ).write_parquet(index_path)
    assets[index_task.asset_id] = FileAsset(index_path)
    protein_tasks = []
    for n, zs in enumerate(z_by_protein):
        task = BuildSlimProteinParquetTask.create(
            index_task=index_task,
            protein=_protein(n),
            asset_id=f"slim_{n}",
            index_name="test",
        )
        path = tmp_path / f"slim_{n}.parquet"
        pl.DataFrame(
            {
                GWASLAB_BETA_COL: [z * _SE for z in zs],
                GWASLAB_SE_COL: [_SE] * len(zs),
            },
            schema={GWASLAB_BETA_COL: pl.Float32, GWASLAB_SE_COL: pl.Float32},
        ).write_parquet(path)
        assets[task.asset_id] = FileAsset(path)
        protein_tasks.append(task)
    task = PppProteinHitTableTask.create(
        asset_id="hits",
        protein_tasks=tuple(protein_tasks),
        index_task=index_task,
        p_threshold=p_threshold,
        merge_gap_bp=merge_gap_bp,
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_significant_variants_within_gap_form_one_segment(tmp_path: Path) -> None:
    variants = [(1, 1_000_000), (1, 1_200_000), (1, 1_400_000)]
    table = _run(tmp_path, variants, [[_HIT_Z, _HIT_Z, _HIT_Z]])
    assert table.height == 1
    row = table.row(0, named=True)
    assert row[HIT_SEG_START_COL] == 1_000_000
    assert row[HIT_SEG_END_COL] == 1_400_000
    assert row[HIT_N_SIG_VARIANTS_COL] == 3
    assert row[HIT_PROTEIN_ID_COL] == "OID0"
    assert row[HIT_GENE_COL] == "G0"
    assert row[HIT_CHR_COL] == 1


def test_gap_larger_than_merge_gap_splits_segments(tmp_path: Path) -> None:
    table = _run(tmp_path, [(1, 1_000_000), (1, 1_300_000)], [[_HIT_Z, _HIT_Z]])
    assert table.sort(HIT_SEG_START_COL)[HIT_SEG_START_COL].to_list() == [
        1_000_000,
        1_300_000,
    ]


def test_segments_do_not_span_chromosomes(tmp_path: Path) -> None:
    table = _run(tmp_path, [(1, 1_000_000), (2, 1_000_010)], [[_HIT_Z, _HIT_Z]])
    assert sorted(table[HIT_CHR_COL].to_list()) == [1, 2]


def test_lead_is_the_min_p_variant_of_the_segment(tmp_path: Path) -> None:
    variants = [(1, 1_000_000), (1, 1_100_000), (1, 1_200_000)]
    table = _run(tmp_path, variants, [[_HIT_Z, 9.0, -7.0]])
    assert table[HIT_LEAD_POS_COL].to_list() == [1_100_000]


def test_nonsignificant_variants_do_not_extend_segments(tmp_path: Path) -> None:
    # The middle null variant would bridge the two hits if it counted.
    variants = [(1, 1_000_000), (1, 1_200_000), (1, 1_400_000)]
    table = _run(tmp_path, variants, [[_HIT_Z, _NULL_Z, _HIT_Z]])
    assert table.height == 2
    assert table[HIT_N_SIG_VARIANTS_COL].to_list() == [1, 1]


def test_threshold_and_gap_are_configurable(tmp_path: Path) -> None:
    variants = [(1, 1_000_000), (1, 1_100_000), (1, 5_000_000)]
    zs = [[_HIT_Z, _HIT_Z, _WEAK_Z]]
    strict = _run(tmp_path / "a", variants, zs)
    lenient = _run(tmp_path / "b", variants, zs, p_threshold=1e-5)
    no_merge = _run(tmp_path / "c", variants, zs, merge_gap_bp=0)
    assert strict.height == 1
    assert lenient.height == 2
    assert no_merge.height == 2


def test_extreme_z_does_not_underflow(tmp_path: Path) -> None:
    table = _run(tmp_path, [(1, 1_000_000)], [[40.0]])
    mlog10p = table[HIT_LEAD_MLOG10P_COL].item()
    assert 300 < mlog10p < float("inf")


def test_protein_without_hits_contributes_no_rows(tmp_path: Path) -> None:
    table = _run(
        tmp_path, [(1, 1_000_000), (1, 2_000_000)], [[_HIT_Z, 0.0], [0.0, 1.0]]
    )
    assert table[HIT_PROTEIN_ID_COL].to_list() == ["OID0"]
