import json
from pathlib import Path, PurePath

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    DIAGNOSTICS_JSON_FILENAME,
    GAMMA_RAW_COL,
    GAMMA_STANDARDIZED_COL,
    WEIGHTS_PARQUET_FILENAME,
    RidgeAnnotationWeightsTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf

# Output of the task on the noisy fixture in test_noisy_shrunk_fit_is_numerically_pinned,
# recorded from a known-good run.
PINNED_GAMMA_RAW = [2.0197461076516072, -0.7935015322177634, 0.656901484778738]
PINNED_GAMMA_STD = [2.026872633214789, -0.7657972111684178, 0.6751477573179109]
PINNED_ALPHA = 10.0
PINNED_INTERCEPT = 2.7887066141341497
PINNED_MEAN_R2 = 0.3693738589448255


def _make_inputs(tmp_path: Path) -> tuple[Path, Path, dict[str, float]]:
    # Use real baseline-LF annotation names (rather than placeholders like
    # "annotA") so family_for_annotation resolves a family for each column
    # instead of raising on an unrecognized name. The annotation matrix and
    # snpvar_meta both carry alleles (A1/A2); the join is allele-aware.
    rng = np.random.default_rng(0)
    n_per_chrom = 400
    truth = {
        "Coding_UCSC_common": 2.0,
        "Promoter_UCSC_common": -1.0,
        "H3K27ac_Hnisz_common": 0.5,
    }

    def _row(chrom: int, bp: int, snp: str, a2: str) -> dict:
        a, b, c = rng.normal(), rng.normal(), rng.normal()
        y = (
            3.0
            + truth["Coding_UCSC_common"] * a
            + truth["Promoter_UCSC_common"] * b
            + truth["H3K27ac_Hnisz_common"] * c
        )
        return {
            "CHR": chrom,
            "BP": bp,
            "SNP": snp,
            "A1": "A",
            "A2": a2,
            "Coding_UCSC_common": a,
            "Promoter_UCSC_common": b,
            "H3K27ac_Hnisz_common": c,
            "snpvar_bin": y,
        }

    rows: list[dict] = []
    for chrom in (1, 2):
        for i in range(n_per_chrom):
            rows.append(_row(chrom, i + 1, f"rs{chrom}_{i}", a2="G"))
    # A multiallelic site: a SECOND allele at chr1 BP 1, sharing the rsID of the
    # first (rs1_0) but with a distinct allele key (A/C vs A/G) and its own
    # annotations + snpvar. The old SNP-keyed+deduped join would have collapsed
    # this to one row; the allele-aware join keeps both, so n_variants == 801.
    rows.append(_row(1, 1, "rs1_0", a2="C"))

    frame = pl.DataFrame(rows)
    annot_path = tmp_path / "annot.parquet"
    frame.drop("snpvar_bin").write_parquet(annot_path)
    meta_path = tmp_path / "meta.parquet"
    frame.select("CHR", "BP", "SNP", "A1", "A2", "snpvar_bin").write_parquet(meta_path)
    return annot_path, meta_path, truth


def _run_task(
    tmp_path: Path, annot_path: Path, meta_path: Path, alphas: tuple[float, ...]
) -> Path:
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    annot_task = FakeTask(
        ReferenceFileMeta(
            group="polyfun",
            sub_group="annotations",
            sub_folder=PurePath("raw"),
            id=AssetId("annot"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
    )
    meta_task = FakeTask(
        SimpleFileMeta("meta", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )
    task = RidgeAnnotationWeightsTask.create(
        asset_id="ridge_weights",
        annotation_parquet_task=annot_task,
        snpvar_meta_task=meta_task,
        alphas=alphas,
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "annot":
            return FileAsset(annot_path)
        if asset_id == "meta":
            return FileAsset(meta_path)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    return result.path


def test_recovers_known_linear_weights(tmp_path: Path):
    annot_path, meta_path, truth = _make_inputs(tmp_path)
    out = _run_task(tmp_path, annot_path, meta_path, alphas=(1e-4, 1e-2, 1.0, 100.0))

    weights = pl.read_parquet(out / WEIGHTS_PARQUET_FILENAME)
    got = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
    for name, expected in truth.items():
        assert abs(got[name] - expected) < 1e-2

    diagnostics = json.loads((out / DIAGNOSTICS_JSON_FILENAME).read_text())
    assert diagnostics["mean_heldout_r2"] > 0.999
    # 800 biallelic + 1 second allele at the multiallelic site, all retained.
    assert diagnostics["n_variants"] == 801


def test_noisy_shrunk_fit_is_numerically_pinned(tmp_path: Path):
    # Regression guard for the ridge arithmetic: a noisy target at a penalty large
    # enough to visibly shrink the coefficients, so the pinned values depend on the
    # standardization, alpha selection, and held-out R^2, not just on the truth.
    annot_path, meta_path, _ = _make_inputs(tmp_path)
    rng = np.random.default_rng(1)
    meta = pl.read_parquet(meta_path)
    meta.with_columns(
        pl.col("snpvar_bin") + pl.Series(rng.normal(scale=3.0, size=meta.height))
    ).write_parquet(meta_path)
    out = _run_task(tmp_path, annot_path, meta_path, alphas=(10.0, 100.0, 1000.0))

    weights = pl.read_parquet(out / WEIGHTS_PARQUET_FILENAME)
    diagnostics = json.loads((out / DIAGNOSTICS_JSON_FILENAME).read_text())
    assert np.allclose(weights[GAMMA_RAW_COL].to_numpy(), PINNED_GAMMA_RAW, atol=1e-10)
    assert np.allclose(
        weights[GAMMA_STANDARDIZED_COL].to_numpy(), PINNED_GAMMA_STD, atol=1e-10
    )
    assert diagnostics["alpha"] == PINNED_ALPHA
    assert diagnostics["intercept"] == pytest.approx(PINNED_INTERCEPT, abs=1e-10)
    assert diagnostics["mean_heldout_r2"] == pytest.approx(PINNED_MEAN_R2, abs=1e-10)
