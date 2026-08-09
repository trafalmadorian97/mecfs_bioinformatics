from pathlib import Path, PurePath

from mecfs_bio.build_system.wf.remote_executor.fake_remote_executor import (
    FakeRemoteExecutor,
)
from mecfs_bio.build_system.wf.remote_executor.remote_job import (
    RemoteJob,
    RemoteResources,
)


def test_fake_records_job_and_writes_declared_outputs(tmp_path: Path) -> None:
    ex = FakeRemoteExecutor(stub_outputs={PurePath("out/pip.txt"): "SNP\tPIP\n"})
    job = RemoteJob(
        image="i",
        commands=["gctb"],
        input_files={},
        s3_inputs={},
        output_files=[PurePath("out/pip.txt")],
        resources=RemoteResources(1, 1, 1),
    )
    ex.run(job, tmp_path)
    assert ex.last_job is job
    assert (tmp_path / "out/pip.txt").read_text().startswith("SNP")
