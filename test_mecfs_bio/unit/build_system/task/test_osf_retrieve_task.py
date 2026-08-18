import hashlib
import shlex
from functools import partial
from pathlib import Path, PurePath
from subprocess import CalledProcessError

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.task.osf_retrieve_task import OSFRetrievalTask
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.util.subproc.run_command import execute_command_with_retries

PAYLOAD = b"chr pos beta\n1 100 0.5\n"
# Passed to the task so a successful download also has to clear verify_hash.
PAYLOAD_MD5 = hashlib.md5(PAYLOAD).hexdigest()


def make_fake_osf_cli(attempts: list[list[str]], failures: int):
    """
    Build a stand-in for the osf CLI that reproduces the two behaviours the
    retry logic has to cope with, both read off osfclient 0.0.5:

    - osfclient opens the local file for writing before it starts downloading,
      so a download that fails part way through still leaves a file behind.
    - without --force, osf fetch refuses to overwrite an existing local file
      and exits non-zero.

    The first `failures` attempts fail the way a transient HTTP 403 from OSF's
    file storage does; the rest succeed.
    """

    def fake_osf_cli(cmd: list[str]) -> str:
        # execute_command joins the command and runs it through a shell, so the
        # fake has to undo the same quoting a shell would.
        args = shlex.split(" ".join(cmd))
        local_path = Path(args[-1])
        if local_path.exists() and "--force" not in args:
            raise CalledProcessError(
                returncode=1,
                cmd=cmd,
                output=f"Local file {local_path} already exists, not overwriting.",
            )
        attempts.append(args)
        local_path.write_bytes(b"")
        if len(attempts) <= failures:
            raise CalledProcessError(
                returncode=1, cmd=cmd, output="Response has status code 403."
            )
        local_path.write_bytes(PAYLOAD)
        return ""

    return fake_osf_cli


def fake_fetch(asset_id: AssetId) -> Asset:
    raise NotImplementedError("OSFRetrievalTask has no dependencies to fetch")


def test_transient_download_failure_is_retried(tmp_path: Path):
    """
    A fetch that fails once and then succeeds produces the asset.

    This exercises the retry as a whole: it only passes if the retried command
    can overwrite the empty file the failed attempt left behind.
    """
    attempts: list[list[str]] = []
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True, exist_ok=True)
    task = OSFRetrievalTask(
        meta=GWASSummaryDataFileMeta(
            id=AssetId("fake_osf_asset"),
            trait="ME_CFS",
            project="DecodeME",
            sub_dir="raw",
            project_path=PurePath("Summary Statistics") / "gwas.regenie.gz",
        ),
        osf_project_id="fake_project",
        md5_hash=PAYLOAD_MD5,
        run_command=partial(
            execute_command_with_retries,
            base_backoff_seconds=0.0,
            sleep=lambda _: None,
            executor=make_fake_osf_cli(attempts, failures=1),
        ),
    )

    asset = task.execute(scratch_dir=scratch, fetch=fake_fetch, wf=make_wf())

    assert asset.path.read_bytes() == PAYLOAD
    assert len(attempts) == 2
