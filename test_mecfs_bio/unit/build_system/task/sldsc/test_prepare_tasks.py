from pathlib import Path

import numpy as np

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sldsc.prepare_ldcts_annotation_set_task import (
    PrepareLdctsAnnotationSetTask,
)
from mecfs_bio.build_system.task.sldsc.prepare_sldsc_baseline_task import (
    PrepareSLDSCBaselineTask,
)
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    read_prepared_annotation_set,
    read_prepared_baseline,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from test_mecfs_bio.unit.build_system.task.sldsc.sldsc_fixture import (
    FIXTURE_BASELINE_PATTERN,
    FIXTURE_LABELS,
    FIXTURE_LDCTS_FILENAME,
    FIXTURE_WEIGHTS_PATTERN,
    write_sldsc_fixture,
)


def _dir_task(asset_id: str) -> FakeTask:
    return FakeTask(
        ReferenceDataDirectoryMeta(
            group="g", sub_group="s", sub_folder=Path("extracted"), id=AssetId(asset_id)
        )
    )


def test_prepare_annotation_set_stores_shared_control_once(tmp_path: Path) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    task = PrepareLdctsAnnotationSetTask.create(
        asset_id="prepared_set",
        source_directory_task=_dir_task("cts"),
        ldcts_filename=FIXTURE_LDCTS_FILENAME,
    )

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == "cts"
        return DirectoryAsset(fixture.cts_dir)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    prepared = read_prepared_annotation_set(result.path)
    assert [m.label for m in prepared.models] == list(FIXTURE_LABELS)
    controls = {m.control_columns for m in prepared.models}
    assert len(controls) == 1  # one shared control set, referenced by every model
    assert (
        len(next(iter(controls))) == 2
    )  # every column of the control prefix is a control
    assert prepared.ld.matrix.shape == (3_000, len(FIXTURE_LABELS) + 2)
    np.testing.assert_allclose(
        prepared.ld.m_5_50.sum(), 2000.0 + 1500.0 + 300.0 + 400.0 + 500.0
    )


def test_prepare_baseline_drops_zero_variance_column_and_unweighted_snps(
    tmp_path: Path,
) -> None:
    fixture = write_sldsc_fixture(tmp_path)
    task = PrepareSLDSCBaselineTask.create(
        asset_id="prepared_baseline",
        baseline_directory_task=_dir_task("baseline"),
        baseline_pattern=FIXTURE_BASELINE_PATTERN,
        weights_directory_task=_dir_task("weights"),
        weights_pattern=FIXTURE_WEIGHTS_PATTERN,
    )

    def fetch(asset_id: AssetId) -> Asset:
        return DirectoryAsset(
            {"baseline": fixture.baseline_dir, "weights": fixture.weights_dir}[asset_id]
        )

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    prepared = read_prepared_baseline(result.path)
    assert not any("ZeroL2" in s for s in prepared.ld.sources)
    assert prepared.ld.matrix.shape == (
        2_970,
        3,
    )  # 30 SNPs without a weight LD score dropped
    assert prepared.w_ld.height == 2_970
