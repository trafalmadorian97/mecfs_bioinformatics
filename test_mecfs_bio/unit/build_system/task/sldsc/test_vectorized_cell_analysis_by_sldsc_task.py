from pathlib import Path
from typing import Protocol

import gwaslab as gl
import numpy as np
import pandas as pd
import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    CellAnalysisByLDSCTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PrepareLdctsAnnotationSetTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_sldsc_baseline_task import (
    PrepareSLDSCBaselineTask,
)
from mecfs_bio.build_system.task.sldsc.vectorized_cell_analysis_by_sldsc_task import (
    VectorizedCellAnalysisBySLDSCTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.sldsc_constants import (
    SLDSC_CTS_COEF_COL,
    SLDSC_CTS_COEF_SE_COL,
    SLDSC_CTS_NAME_COL,
    SLDSC_CTS_P_COL,
    SLDSC_DEFAULT_MIN_HAPMAP3_SNPS,
    SLDSCDuplicateRsidPolicy,
)
from test_mecfs_bio.unit.build_system.task.sldsc.sldsc_fixture import (
    FIXTURE_BASELINE_PATTERN,
    FIXTURE_LABELS,
    FIXTURE_LDCTS_FILENAME,
    FIXTURE_WEIGHTS_PATTERN,
    SLDSCFixture,
    sumstats_from_frame,
    write_sldsc_fixture,
)


def _dir_task(asset_id: str) -> FakeTask:
    return FakeTask(
        ReferenceDataDirectoryMeta(
            group="g", sub_group="s", sub_folder=Path("extracted"), id=AssetId(asset_id)
        )
    )


_SUMSTATS = FakeTask(
    GWASLabSumStatsMeta(id=AssetId("sumstats"), trait="t", project="p")
)
_BASELINE_DIR = _dir_task("baseline_dir")
_WEIGHTS_DIR = _dir_task("weights_dir")
_CTS_DIR = _dir_task("cts_dir")
_PREP_BASELINE = PrepareSLDSCBaselineTask.create(
    asset_id="prep_baseline",
    baseline_directory_task=_BASELINE_DIR,
    baseline_pattern=FIXTURE_BASELINE_PATTERN,
    weights_directory_task=_WEIGHTS_DIR,
    weights_pattern=FIXTURE_WEIGHTS_PATTERN,
)
_PREP_SET = PrepareLdctsAnnotationSetTask.create(
    asset_id="prep_set",
    source_directory_task=_CTS_DIR,
    ldcts_filename=FIXTURE_LDCTS_FILENAME,
)


@frozen(slots=True)
class _Store:
    """Materializes fixture assets and runs Tasks against them, like the build system's fetch."""

    root: Path
    assets: dict[str, Asset]

    def fetch(self, asset_id: AssetId) -> Asset:
        return self.assets[str(asset_id)]

    def run(self, task: Task) -> Asset:
        scratch = self.root / f"scratch_{task.asset_id}"
        scratch.mkdir()
        asset = task.execute(scratch_dir=scratch, fetch=self.fetch, wf=make_wf())
        self.assets[str(task.asset_id)] = asset
        return asset


def _store(
    tmp_path: Path, fixture: SLDSCFixture, sumstats_frame: pd.DataFrame
) -> _Store:
    tmp_path.mkdir(parents=True, exist_ok=True)
    pickle_path = tmp_path / "sumstats.pickle"
    gl.dump_pickle(sumstats_from_frame(sumstats_frame), path=str(pickle_path))
    store = _Store(
        root=tmp_path,
        assets={
            "sumstats": FileAsset(pickle_path),
            "baseline_dir": DirectoryAsset(fixture.baseline_dir),
            "weights_dir": DirectoryAsset(fixture.weights_dir),
            "cts_dir": DirectoryAsset(fixture.cts_dir),
        },
    )
    store.run(_PREP_BASELINE)
    store.run(_PREP_SET)
    return store


# The fixture has ~2,950 HapMap3-style regression SNPs, far below the production default.
_TEST_MIN_HAPMAP3_SNPS = 1_000


def _vectorized(
    asset_id: str = "vec",
    min_hapmap3_snps: int = _TEST_MIN_HAPMAP3_SNPS,
    duplicate_rsids: SLDSCDuplicateRsidPolicy = "raise",
) -> VectorizedCellAnalysisBySLDSCTask:
    return VectorizedCellAnalysisBySLDSCTask.create(
        asset_id=asset_id,
        source_sumstats_task=_SUMSTATS,
        baseline_task=_PREP_BASELINE,
        annotation_set_task=_PREP_SET,
        min_hapmap3_snps=min_hapmap3_snps,
        duplicate_rsids=duplicate_rsids,
    )


def _read_result(asset: Asset) -> pl.DataFrame:
    assert isinstance(asset, FileAsset)
    return pl.read_csv(asset.path).sort(SLDSC_CTS_NAME_COL)


def test_exact_mode_matches_gwaslab_cell_type_specific(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, fixture.sumstats_frame)
    old = CellAnalysisByLDSCTask.create(
        asset_id="old",
        source_sumstats_task=_SUMSTATS,
        ref_ld_chr_cts_task=_CTS_DIR,
        ref_ld_chr_cts_filename=FIXTURE_LDCTS_FILENAME,
        ref_ld_chr_task=_BASELINE_DIR,
        ref_ld_chr_inner_dirname=FIXTURE_BASELINE_PATTERN,
        w_ld_chr_task=_WEIGHTS_DIR,
        w_ld_chr_inner_dirname=FIXTURE_WEIGHTS_PATTERN,
    )
    expected = _read_result(store.run(old))
    actual = _read_result(store.run(_vectorized()))
    assert actual[SLDSC_CTS_NAME_COL].to_list() == sorted(FIXTURE_LABELS)
    for col in (SLDSC_CTS_COEF_COL, SLDSC_CTS_COEF_SE_COL, SLDSC_CTS_P_COL):
        np.testing.assert_allclose(
            actual[col].to_numpy(), expected[col].to_numpy(), rtol=1e-5
        )


def test_degenerate_z_rows_do_not_change_results(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    clean = _read_result(
        _store(tmp_path / "clean", fixture, fixture.sumstats_frame).run(_vectorized())
    )
    degenerate = fixture.sumstats_frame.copy()
    degenerate.loc[[10, 11], "SE"] = 0.0  # +-inf Z: removed by the chi-square filter
    degenerate.loc[[12], ["BETA", "SE"]] = 0.0  # NaN Z: dropped
    with_bad = _read_result(
        _store(tmp_path / "degenerate", fixture, degenerate).run(_vectorized())
    )
    without_rows = fixture.sumstats_frame.drop(index=[10, 11, 12])
    dropped = _read_result(
        _store(tmp_path / "dropped", fixture, without_rows).run(_vectorized())
    )
    np.testing.assert_allclose(
        with_bad[SLDSC_CTS_COEF_COL].to_numpy(),
        dropped[SLDSC_CTS_COEF_COL].to_numpy(),
        rtol=1e-12,
    )
    assert clean.height == with_bad.height


class _FrameCorruption(Protocol):
    def __call__(self, frame: pd.DataFrame) -> pd.DataFrame: ...


def _missing_n_on_regression_snp(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.assign(N=frame["N"].astype(float).where(frame.index != 500))


def _duplicate_rsid(frame: pd.DataFrame) -> pd.DataFrame:
    return pd.concat([frame, frame.iloc[[400]]], ignore_index=True)


def _no_n_column(frame: pd.DataFrame) -> pd.DataFrame:
    return frame.drop(columns="N")


@pytest.mark.parametrize(
    "corrupt", [_missing_n_on_regression_snp, _duplicate_rsid, _no_n_column]
)
def test_invalid_sumstats_fail_fast(tmp_path: Path, corrupt: _FrameCorruption) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, corrupt(fixture.sumstats_frame))
    with pytest.raises(AssertionError):
        store.run(_vectorized())


def _assert_same_coefficients(actual: pl.DataFrame, expected: pl.DataFrame) -> None:
    for col in (SLDSC_CTS_COEF_COL, SLDSC_CTS_COEF_SE_COL):
        np.testing.assert_allclose(
            actual[col].to_numpy(), expected[col].to_numpy(), rtol=1e-12
        )


def test_drop_duplicate_rsids_drops_every_copy(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    frame = fixture.sumstats_frame
    conflicting = frame.iloc[[400]].assign(BETA=-frame["BETA"].iloc[400])
    duplicated = pd.concat([frame, conflicting], ignore_index=True)
    actual = _read_result(
        _store(tmp_path / "duplicated", fixture, duplicated).run(
            _vectorized(duplicate_rsids="drop")
        )
    )
    expected = _read_result(
        _store(tmp_path / "removed", fixture, frame.drop(index=[400])).run(
            _vectorized()
        )
    )
    _assert_same_coefficients(actual, expected)


def test_drop_duplicate_rsids_keeps_snp_whose_duplicate_fails_allele_filter(
    tmp_path: Path,
) -> None:
    """
    A row sharing a regression SNP's rsID but carrying alleles that do not match HapMap3 (for
    example a multi-allelic sibling) is removed by the HapMap3 allele filter, so the regression
    SNP itself is not a duplicate and is kept.
    """
    fixture = write_sldsc_fixture(tmp_path)
    frame = fixture.sumstats_frame
    alleles = {frame["EA"].iloc[400], frame["NEA"].iloc[400]}
    other = next(b for b in "ACGT" if b not in alleles)
    sibling = frame.iloc[[400]].assign(EA=other)
    with_sibling = pd.concat([frame, sibling], ignore_index=True)
    actual = _read_result(
        _store(tmp_path / "sibling", fixture, with_sibling).run(
            _vectorized(duplicate_rsids="drop")
        )
    )
    expected = _read_result(
        _store(tmp_path / "clean", fixture, frame).run(_vectorized())
    )
    _assert_same_coefficients(actual, expected)


def test_annotation_set_missing_regression_snps_fails(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, fixture.sumstats_frame)
    small_cts = tmp_path / "small_cts"
    write_sldsc_fixture(
        small_cts
    )  # same SNP ids, then drop chromosome 3 from the cts files
    for f in (small_cts / "cts_root" / "cts").glob("*.3.l2.*"):
        f.unlink()
    store.assets["small_cts_dir"] = DirectoryAsset(small_cts / "cts_root")
    small_set = PrepareLdctsAnnotationSetTask.create(
        asset_id="small_set",
        source_directory_task=_dir_task("small_cts_dir"),
        ldcts_filename=FIXTURE_LDCTS_FILENAME,
    )
    store.run(small_set)
    task = VectorizedCellAnalysisBySLDSCTask.create(
        asset_id="vec_small",
        source_sumstats_task=_SUMSTATS,
        baseline_task=_PREP_BASELINE,
        annotation_set_task=small_set,
        min_hapmap3_snps=_TEST_MIN_HAPMAP3_SNPS,
    )
    with pytest.raises(AssertionError):
        store.run(task)


def test_too_few_hapmap3_snps_fails(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    store = _store(tmp_path, fixture, fixture.sumstats_frame)
    with pytest.raises(AssertionError):
        store.run(_vectorized(min_hapmap3_snps=SLDSC_DEFAULT_MIN_HAPMAP3_SNPS))
