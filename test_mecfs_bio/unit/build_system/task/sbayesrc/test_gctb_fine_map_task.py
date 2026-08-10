import json
from pathlib import Path, PurePath

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.sbayesrc import gctb_gwfm_constants as c
from mecfs_bio.build_system.task.sbayesrc.gctb_fine_map_task import GctbFineMapTask
from mecfs_bio.build_system.task.sbayesrc.gctb_mcmc_options import GctbMcmcOptions
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.build_system.wf.remote_executor.fake_remote_executor import (
    FakeRemoteExecutor,
)


def _dispatch_and_capture_commands(
    tmp_path: Path, mcmc_options: GctbMcmcOptions
) -> list[str]:
    """Build and execute a task with the given MCMC options, returning its commands."""
    ma_path = tmp_path / "trait.ma"
    ma_path.write_text("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")
    ma_meta = FilteredGWASDataMeta(
        id=AssetId("ma_src"),
        trait="mecfs",
        project="decodeme",
        sub_dir="x",
        extension=".ma",
    )
    ma_stub = FakeTask(meta=ma_meta)
    marker_path = tmp_path / "marker.json"
    marker_path.write_text(json.dumps({c.MARKER_PREFIX_KEY: "s3://bucket/ref/v1"}))
    marker_stub = FakeTask(meta=SimpleFileMeta(AssetId("marker")))

    def _fetch(asset_id: AssetId) -> Asset:
        if asset_id == ma_stub.asset_id:
            return FileAsset(ma_path)
        assert asset_id == marker_stub.asset_id
        return FileAsset(marker_path)

    fake = FakeRemoteExecutor(
        stub_outputs={PurePath("work/out/gwfm.snpRes"): "SNP\tPIP\n"}
    )
    task = GctbFineMapTask.create(
        id="gwfm_out",
        ma_task=ma_stub,
        reference_task=marker_stub,
        image="gctb:test",
        mcmc_options=mcmc_options,
    )
    task.execute(scratch_dir=tmp_path, fetch=_fetch, wf=make_wf(remote_executor=fake))
    assert fake.last_job is not None
    return list(fake.last_job.commands)


def _gwfm_command(commands: list[str]) -> str:
    (gwfm,) = [cmd for cmd in commands if "--gwfm RC" in cmd]
    return gwfm


def test_mcmc_options_are_appended_to_the_gwfm_step(tmp_path: Path) -> None:
    commands = _dispatch_and_capture_commands(
        tmp_path, GctbMcmcOptions(chain_length=100, burn_in=20, seed=123)
    )
    gwfm = _gwfm_command(commands)
    assert "--chain-length 100" in gwfm
    assert "--burn-in 20" in gwfm
    assert "--seed 123" in gwfm


def test_default_mcmc_options_leave_the_gwfm_step_flagless(tmp_path: Path) -> None:
    commands = _dispatch_and_capture_commands(tmp_path, GctbMcmcOptions())
    gwfm = _gwfm_command(commands)
    assert "--chain-length" not in gwfm
    assert "--burn-in" not in gwfm
    # Only the gwfm step carries MCMC overrides; the other steps never do.
    assert not any("--chain-length" in cmd for cmd in commands)


def test_dispatches_a_wellformed_gwfm_job(tmp_path: Path) -> None:
    ma_path = tmp_path / "trait.ma"
    ma_path.write_text("SNP\tA1\tA2\tfreq\tb\tse\tp\tN\n")
    ma_meta = FilteredGWASDataMeta(
        id=AssetId("ma_src"),
        trait="mecfs",
        project="decodeme",
        sub_dir="x",
        extension=".ma",
    )
    ma_stub = FakeTask(meta=ma_meta)

    marker_path = tmp_path / "marker.json"
    s3_prefix = "s3://bucket/sbayesrc/ld/Imputed13M/v1"
    marker_path.write_text(json.dumps({c.MARKER_PREFIX_KEY: s3_prefix}))
    marker_stub = FakeTask(meta=SimpleFileMeta(AssetId("marker")))

    def _fetch(asset_id: AssetId) -> Asset:
        if asset_id == ma_stub.asset_id:
            return FileAsset(ma_path)
        assert asset_id == marker_stub.asset_id
        return FileAsset(marker_path)

    fake = FakeRemoteExecutor(
        stub_outputs={PurePath("work/out/gwfm.snpRes"): "SNP\tPIP\n"}
    )
    wf = make_wf(remote_executor=fake)

    task = GctbFineMapTask.create(
        id="gwfm_out",
        ma_task=ma_stub,
        reference_task=marker_stub,
        image="gctb:test",
    )
    asset = task.execute(scratch_dir=tmp_path, fetch=_fetch, wf=wf)

    job = fake.last_job
    assert job is not None
    assert any("--gwfm RC" in cmd for cmd in job.commands)
    assert any(str(remote).endswith(".ma") for remote in job.input_files.values())
    assert job.s3_inputs == {s3_prefix: PurePath(c.REMOTE_REF_DIR)}
    assert job.resources.memory_gb == c.DEFAULT_MEMORY_GB
    assert isinstance(asset, DirectoryAsset)


def test_deps_are_the_ma_and_reference_tasks() -> None:
    ma_meta = FilteredGWASDataMeta(
        id=AssetId("ma_src"), trait="t", project="p", sub_dir="x", extension=".ma"
    )
    ma_stub = FakeTask(meta=ma_meta)
    marker_stub = FakeTask(meta=SimpleFileMeta(AssetId("marker")))
    task = GctbFineMapTask.create(
        id="gwfm_out", ma_task=ma_stub, reference_task=marker_stub, image="gctb:test"
    )
    assert task.deps == [ma_stub, marker_stub]
