"""Synthetic inputs shared by the fast CT-LDSC task tests: a HapMap3-like index with LD scores,
and trait sumstats built from known z-scores on that index."""

from pathlib import Path, PurePath

import gwaslab as gl
import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.base_task import Task
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
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    HAPMAP3_SNPLIST_A1_COL,
    HAPMAP3_SNPLIST_A2_COL,
    HAPMAP3_SNPLIST_RSID_COL,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

# Non-palindromic allele pairs, so no variant is dropped as strand-ambiguous.
_ALLELE_PAIRS = (("A", "G"), ("C", "T"), ("A", "C"), ("G", "T"))
TRAIT_SE = 0.1


@frozen(slots=True)
class SyntheticIndex:
    """The index task, the assets every fetch needs, and the index variants (rsID, EA, NEA, CHR,
    POS, L2) in genome order."""

    index_task: FastCtLdscVariantIndexTask
    assets: dict[AssetId, Asset]
    variants: pl.DataFrame


def _reference_meta(asset_id: str) -> ReferenceFileMeta:
    return ReferenceFileMeta(
        group="g",
        sub_group="s",
        sub_folder=PurePath("raw"),
        id=AssetId(asset_id),
        extension=".parquet",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
    )


def build_synthetic_index(tmp_path: Path, n_snps: int, seed: int) -> SyntheticIndex:
    """n_snps variants split over chromosomes 1 and 2 (far from the MHC), with LD scores drawn
    uniformly on [1, 50]. The index task is executed so its asset is in the returned assets."""
    rng = np.random.default_rng(seed)
    half = n_snps // 2
    chrom = np.where(np.arange(n_snps) < half, 1, 2)
    pos = 1_000 * (np.arange(n_snps) % max(half, 1) + 1)
    pairs = [_ALLELE_PAIRS[i % len(_ALLELE_PAIRS)] for i in range(n_snps)]
    variants = pl.DataFrame(
        {
            GWASLAB_RSID_COL: [f"rs{i}" for i in range(n_snps)],
            GWASLAB_EFFECT_ALLELE_COL: [p[0] for p in pairs],
            GWASLAB_NON_EFFECT_ALLELE_COL: [p[1] for p in pairs],
            GWASLAB_CHROM_COL: chrom.astype(np.int64),
            GWASLAB_POS_COL: pos.astype(np.int64),
            LD_SCORE_LD_SCORE_COL: rng.uniform(1.0, 50.0, size=n_snps),
        }
    )
    snplist = variants.select(
        pl.col(GWASLAB_RSID_COL).alias(HAPMAP3_SNPLIST_RSID_COL),
        pl.col(GWASLAB_EFFECT_ALLELE_COL).alias(HAPMAP3_SNPLIST_A1_COL),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(HAPMAP3_SNPLIST_A2_COL),
    )
    ld = variants.select(
        pl.col(GWASLAB_CHROM_COL).alias(LD_SCORE_CHROM_COL),
        pl.col(GWASLAB_RSID_COL).alias(LD_SCORE_RSID_COL),
        pl.col(GWASLAB_POS_COL).alias(LD_SCORE_POS_COL),
        pl.col(LD_SCORE_LD_SCORE_COL),
        pl.lit(float(n_snps)).alias(LD_SCORE_M_5_50_COL),
    )
    hapmap3_task = FakeTask(meta=_reference_meta("hapmap3"))
    ld_task = FakeTask(meta=_reference_meta("ld_scores"))
    assets: dict[AssetId, Asset] = {}
    register_frame(tmp_path, assets, hapmap3_task, snplist)
    register_frame(tmp_path, assets, ld_task, ld)
    index_task = FastCtLdscVariantIndexTask.create(
        "index", hapmap3_task=hapmap3_task, ld_scores_task=ld_task
    )
    assets[index_task.asset_id] = run_task(tmp_path, index_task, assets)
    return SyntheticIndex(index_task=index_task, assets=assets, variants=variants)


def register_frame(
    tmp_path: Path, assets: dict[AssetId, Asset], task: Task, frame: pl.DataFrame
) -> None:
    path = tmp_path / f"{task.asset_id}.parquet"
    frame.write_parquet(path)
    assets[task.asset_id] = FileAsset(path)


def register_pickle(
    tmp_path: Path, assets: dict[AssetId, Asset], task: Task, sumstats: gl.Sumstats
) -> None:
    path = tmp_path / f"{task.asset_id}.pickle"
    gl.dump_pickle(sumstats, path=str(path))
    assets[task.asset_id] = FileAsset(path)


def run_task(tmp_path: Path, task: Task, assets: dict[AssetId, Asset]) -> Asset:
    scratch = tmp_path / f"scratch_{task.asset_id}"
    scratch.mkdir()
    return task.execute(scratch_dir=scratch, fetch=lambda a: assets[a], wf=make_wf())


def trait_frame(
    variants: pl.DataFrame,
    z: np.ndarray,
    *,
    n: np.ndarray | float,
    flip_every: int = 2,
) -> pl.DataFrame:
    """A sumstats frame carrying the given index-oriented z on every index variant. Every
    flip_every-th row reports the swapped allele orientation with BETA negated, so alignment must
    flip it back. NaN entries of z are omitted from the frame."""
    flip = np.arange(variants.height) % flip_every == 0
    beta = np.where(flip, -z, z) * TRAIT_SE
    frame = pl.DataFrame(
        {
            GWASLAB_RSID_COL: variants[GWASLAB_RSID_COL],
            GWASLAB_CHROM_COL: variants[GWASLAB_CHROM_COL],
            GWASLAB_POS_COL: variants[GWASLAB_POS_COL],
            GWASLAB_EFFECT_ALLELE_COL: np.where(
                flip,
                variants[GWASLAB_NON_EFFECT_ALLELE_COL],
                variants[GWASLAB_EFFECT_ALLELE_COL],
            ),
            GWASLAB_NON_EFFECT_ALLELE_COL: np.where(
                flip,
                variants[GWASLAB_EFFECT_ALLELE_COL],
                variants[GWASLAB_NON_EFFECT_ALLELE_COL],
            ),
            GWASLAB_BETA_COL: beta,
            GWASLAB_SE_COL: np.full(variants.height, TRAIT_SE),
            GWASLAB_P_COL: np.full(variants.height, 0.5),
            GWASLAB_SAMPLE_SIZE_COLUMN: np.broadcast_to(n, (variants.height,)),
        }
    )
    return frame.filter(pl.Series(np.isfinite(z)))


def trait_sumstats(frame: pl.DataFrame, build: GenomeBuild | None) -> gl.Sumstats:
    """The frame as a gwaslab Sumstats object (EA/NEA become pandas Categoricals)."""
    return gl.Sumstats(
        frame.to_pandas(),
        rsid=GWASLAB_RSID_COL,
        chrom=GWASLAB_CHROM_COL,
        pos=GWASLAB_POS_COL,
        ea=GWASLAB_EFFECT_ALLELE_COL,
        nea=GWASLAB_NON_EFFECT_ALLELE_COL,
        beta=GWASLAB_BETA_COL,
        se=GWASLAB_SE_COL,
        p=GWASLAB_P_COL,
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
        build=build if build is not None else "99",
        verbose=False,
    )


def parquet_trait_meta(
    asset_id: str, genome_build: GenomeBuild | None
) -> FilteredGWASDataMeta:
    return FilteredGWASDataMeta(
        id=AssetId(asset_id),
        trait=asset_id,
        project="test",
        sub_dir="processed",
        read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        genome_build=genome_build,
    )


def pickle_trait_meta(
    asset_id: str, genome_build: GenomeBuild | None
) -> GWASLabSumStatsMeta:
    return GWASLabSumStatsMeta(
        id=AssetId(asset_id), trait=asset_id, project="test", genome_build=genome_build
    )
