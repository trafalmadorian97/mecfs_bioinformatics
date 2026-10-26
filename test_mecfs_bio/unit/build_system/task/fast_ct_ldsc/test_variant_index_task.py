"""The shared HapMap3 x LD-score variant index and the regression context built from it."""

from pathlib import Path, PurePath

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.consolidate_ld_scores_task import (
    LD_SCORE_CHROM_COL,
    LD_SCORE_LD_SCORE_COL,
    LD_SCORE_M_5_50_COL,
    LD_SCORE_POS_COL,
    LD_SCORE_RSID_COL,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fast_ct_ldsc.variant_index_task import (
    FastCtLdscVariantIndexTask,
    load_indexed_ldsc_context,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    HAPMAP3_SNPLIST_A1_COL,
    HAPMAP3_SNPLIST_A2_COL,
    HAPMAP3_SNPLIST_RSID_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

# Inside the hg19 extended MHC (chr6:25.7-33.4 Mb).
_MHC_POS = 30_000_000


def _meta(asset_id: str) -> ReferenceFileMeta:
    return ReferenceFileMeta(
        group="g",
        sub_group="s",
        sub_folder=PurePath("raw"),
        id=AssetId(asset_id),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def _snplist() -> pl.DataFrame:
    # rs_snplist_only is absent from the LD scores; rs_pal is strand-ambiguous (A/T).
    return pl.DataFrame(
        {
            HAPMAP3_SNPLIST_RSID_COL: [
                "rs3",
                "rs1",
                "rs_mhc",
                "rs_pal",
                "rs_snplist_only",
            ],
            HAPMAP3_SNPLIST_A1_COL: ["C", "A", "G", "A", "C"],
            HAPMAP3_SNPLIST_A2_COL: ["T", "G", "A", "T", "G"],
        }
    )


def _ld() -> pl.DataFrame:
    # rs_ld_only is absent from the snplist.
    return pl.DataFrame(
        {
            LD_SCORE_CHROM_COL: [1, 1, 2, 6, 1],
            LD_SCORE_RSID_COL: ["rs1", "rs_pal", "rs3", "rs_mhc", "rs_ld_only"],
            LD_SCORE_POS_COL: [100, 150, 50, _MHC_POS, 200],
            LD_SCORE_LD_SCORE_COL: [1.0, 2.0, 3.0, 4.0, 5.0],
            LD_SCORE_M_5_50_COL: [600.0, 600.0, 300.0, 100.0, 600.0],
        }
    )


def _setup(tmp_path: Path) -> tuple[FastCtLdscVariantIndexTask, dict[AssetId, Asset]]:
    hapmap3 = FakeTask(meta=_meta("hapmap3"))
    ld = FakeTask(meta=_meta("ld"))
    assets: dict[AssetId, Asset] = {}
    for task, frame in ((hapmap3, _snplist()), (ld, _ld())):
        path = tmp_path / f"{task.asset_id}.parquet"
        frame.write_parquet(path)
        assets[task.asset_id] = FileAsset(path)
    index_task = FastCtLdscVariantIndexTask.create(
        "index", hapmap3_task=hapmap3, ld_scores_task=ld
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    index_asset = index_task.execute(
        scratch_dir=scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assets[index_task.asset_id] = index_asset
    return index_task, assets


def test_index_joins_snplist_to_ld_scores_with_a1_as_effect_allele(
    tmp_path: Path,
) -> None:
    index_task, assets = _setup(tmp_path)
    asset = assets[index_task.asset_id]
    assert isinstance(asset, FileAsset)
    out = pl.read_parquet(asset.path)
    # Sorted by (CHR, POS); only rsIDs in both inputs survive.
    assert out[GWASLAB_RSID_COL].to_list() == ["rs1", "rs_pal", "rs3", "rs_mhc"]
    assert out[GWASLAB_EFFECT_ALLELE_COL].to_list() == ["A", "A", "C", "G"]
    assert out[GWASLAB_NON_EFFECT_ALLELE_COL].to_list() == ["G", "T", "T", "A"]
    assert out[GWASLAB_CHROM_COL].to_list() == [1, 1, 2, 6]
    assert out[GWASLAB_POS_COL].to_list() == [100, 150, 50, _MHC_POS]


def test_loaded_context_drops_hg19_mhc_and_strand_ambiguous(tmp_path: Path) -> None:
    index_task, assets = _setup(tmp_path)
    ictx = load_indexed_ldsc_context(index_task, lambda asset_id: assets[asset_id])
    assert ictx.variants[GWASLAB_RSID_COL].to_list() == ["rs1", "rs3"]
    assert ictx.context.ld.tolist() == [1.0, 3.0]
    # M counts every reference variant, before the join.
    assert ictx.context.m == 1000.0
