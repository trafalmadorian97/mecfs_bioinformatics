"""
End-to-end GWFM system test: drives GctbFineMapTask.execute inside the self-built gctb Docker
image via LocalDockerRemoteExecutor on a tiny known-truth toy, asserting the three planted causal
SNPs are recovered with high posterior inclusion probability. The toy reference, sumstats,
annotation and gene map are committed under test_data/gwfm_toy; see
experiments/claude/design_specs/2026-08-09-gwfm-toy-recipe.md for how they were generated.
"""

import json
from pathlib import Path

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c
from mecfs_bio.build_system.task.sbayesrc.gctb_fine_map_task import GctbFineMapTask
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.remote_executor.local_docker_remote_executor import (
    LocalDockerRemoteExecutor,
)
from mecfs_bio.util.subproc.run_command import execute_command

_IMAGE_TAG = "gctb:test"
_TOY_DATA = Path("test_mecfs_bio/system/test_data/gwfm_toy")
_CAUSAL_SNPS = ("rs11", "rs171", "rs331")
_PIP_THRESHOLD = 0.9


def _read_pip(snp_res: Path) -> dict[str, float]:
    lines = snp_res.read_text().splitlines()
    header = lines[0].split()
    name_i, pip_i = header.index("Name"), header.index("PIP")
    return {
        r[name_i]: float(r[pip_i]) for r in (line.split() for line in lines[1:]) if r
    }


def test_gwfm_recovers_planted_causal_snps(tmp_path: Path) -> None:
    execute_command(
        [
            "docker",
            "build",
            "--build-arg",
            f"GCTB_URL={c.GCTB_BINARY_URL}",
            "--build-arg",
            f"GCTB_SHA256={c.GCTB_BINARY_SHA256}",
            "-t",
            _IMAGE_TAG,
            "docker/gctb",
        ]
    )

    ma_path = (_TOY_DATA / "toy.ma").resolve()
    reference_dir = (_TOY_DATA / "reference").resolve()
    marker_path = tmp_path / "marker.json"
    marker_path.write_text(json.dumps({c.MARKER_PREFIX_KEY: str(reference_dir)}))

    ma_meta = FilteredGWASDataMeta(
        id=AssetId("toy_ma"),
        trait="toy",
        project="gwfm",
        sub_dir="gwfm",
        extension=".ma",
        read_spec=None,
    )
    ma_task = FakeTask(meta=ma_meta)
    reference_task = FakeTask.create_with_filemeta("gwfm_reference_marker")

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == ma_task.asset_id:
            return FileAsset(ma_path)
        assert asset_id == reference_task.asset_id
        return FileAsset(marker_path)

    task = GctbFineMapTask.create(
        id="toy_gwfm",
        ma_task=ma_task,
        reference_task=reference_task,
        image=_IMAGE_TAG,
        threads=1,
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    asset = task.execute(
        scratch_dir=scratch,
        fetch=fetch,
        wf=make_wf(remote_executor=LocalDockerRemoteExecutor()),
    )

    assert isinstance(asset, DirectoryAsset)
    pip = _read_pip(asset.path / "gwfm.snpRes")
    for snp in _CAUSAL_SNPS:
        assert pip[snp] > _PIP_THRESHOLD
