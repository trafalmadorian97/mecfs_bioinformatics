import shutil
from pathlib import Path

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_table_task import (
    PolyfunPriorComparisonTableTask,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_variants import (
    CS_PF_EXT_COL,
    CS_U_COL,
    LIFT_EXT_COL,
    LIFT_INT_COL,
    PIP_PF_INT_COL,
    PriorComparisonRuns,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    FILTERED_GWAS_FILENAME,
    PIP_FILENAME,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.polyfun_explain_display_columns import (
    DISP_POS,
)
from test_mecfs_bio.unit.build_system.task.polyfun_explain.test_polyfun_explain_contrast_task import (
    _ExplainInputs,
    build_synthetic_explain_inputs,
)


def _execute_table(
    tmp_path: Path, uni_task: Task, inputs: _ExplainInputs, assets: dict[str, Asset]
) -> Asset:
    # The fixture has one prior run; using it as both the external and the
    # internal prior makes the two prior column sets identical and easy to check.
    task = PolyfunPriorComparisonTableTask.create(
        asset_id="table",
        runs=PriorComparisonRuns(
            susie_uniform_task=uni_task,
            external_prior_susie_task=inputs.pf_task,
            external_prior_contrast_task=inputs.contrast_task,
            internal_prior_susie_task=inputs.pf_task,
            internal_prior_contrast_task=inputs.contrast_task,
        ),
        secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
    )

    def fetch(asset_id: AssetId) -> Asset:
        return assets[str(asset_id)]

    scratch = tmp_path / "table_scratch"
    scratch.mkdir()
    return task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())


def test_table_rows_are_credible_set_union_with_lifts(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)
    result = _execute_table(
        tmp_path, inputs.uni_task, inputs, assets=dict(inputs.fetch_map)
    )
    assert isinstance(result, FileAsset)
    table = pl.read_parquet(result.path)

    # Uniform credible set: POS 10-40; polyfun credible set: POS 10 only. The
    # polyfun run's top variant (POS 10) sorts first.
    assert sorted(table[DISP_POS].to_list()) == [10, 20, 30, 40]
    top = table.row(0, named=True)
    assert top[DISP_POS] == 10
    # The fixture SNPID position is POS + 1000.
    assert top["pos_hg38"] == 1010
    assert top[CS_PF_EXT_COL] == 1
    assert top[CS_U_COL] == 1
    # lift = m * w_i / sum(w) with weights (8, 1, 1, 1, 1, 1).
    assert top[LIFT_EXT_COL] == pytest.approx(6 * 8 / 13)
    assert top[LIFT_INT_COL] == pytest.approx(6 * 8 / 13)
    assert top[PIP_PF_INT_COL] == pytest.approx(0.8)
    # A uniform-only credible-set variant has no polyfun credible-set number.
    assert table.filter(pl.col(DISP_POS) == 20)[CS_PF_EXT_COL].to_list() == [None]


def test_prior_run_with_a_variant_the_uniform_run_lacks_fails(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)
    assets = dict(inputs.fetch_map)
    # Copy the uniform run without its last variant (POS 60, in no credible set),
    # so the prior run carries one variant the uniform run does not.
    uni_dir = assets[str(inputs.uni_task.asset_id)]
    assert isinstance(uni_dir, DirectoryAsset)
    short_dir = tmp_path / "short_uniform"
    shutil.copytree(uni_dir.path, short_dir)
    for filename in (FILTERED_GWAS_FILENAME, PIP_FILENAME):
        path = short_dir / filename
        frame = pl.read_parquet(path)
        frame.head(frame.height - 1).write_parquet(path)
    short_task = FakeTask(
        ResultDirectoryMeta(id=AssetId("short_uniform"), trait="t", project="p")
    )
    assets["short_uniform"] = DirectoryAsset(short_dir)
    with pytest.raises(AssertionError):
        _execute_table(tmp_path, short_task, inputs, assets=assets)
