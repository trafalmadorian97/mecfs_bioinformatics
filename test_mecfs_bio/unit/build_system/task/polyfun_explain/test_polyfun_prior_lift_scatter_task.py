from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_variants import (
    PriorComparisonRuns,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_lift_scatter_task import (
    PolyfunPriorLiftScatterTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from test_mecfs_bio.unit.build_system.task.polyfun_explain.test_polyfun_explain_contrast_task import (
    build_synthetic_explain_inputs,
)


def test_scatter_writes_html(tmp_path: Path):
    # The fixture has one prior run; it serves as both the external and the
    # internal prior.
    inputs = build_synthetic_explain_inputs(tmp_path)
    task = PolyfunPriorLiftScatterTask.create(
        asset_id="scatter",
        runs=PriorComparisonRuns(
            susie_uniform_task=inputs.uni_task,
            external_prior_susie_task=inputs.pf_task,
            external_prior_contrast_task=inputs.contrast_task,
            internal_prior_susie_task=inputs.pf_task,
            internal_prior_contrast_task=inputs.contrast_task,
        ),
        secondary_position=SecondaryPositionFromSnpid(build_label="hg38"),
    )
    fetch_map = dict(inputs.fetch_map)

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    scratch = tmp_path / "scatter_scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    assert result.path.suffix == ".html"
    assert result.path.stat().st_size > 0
