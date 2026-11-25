from pathlib import Path

import attrs
import pytest

from mecfs_bio.util.download.robust_download import (
    ARIA2_MAX_CONNECTIONS_PER_SERVER,
    AriaDownloader,
    DownloadAttemptFailedError,
    Downloader,
    robust_download,
)
from mecfs_bio.util.download.verify import calc_md5_checksum

_TRUE_FILE_CONTENTS = "ABC"
_TRUE_FILE_MD5 = "902fbdd2b1df0c4f70b4a5d23525e932"


def _no_sleep(seconds: float) -> None:
    pass


@attrs.define(slots=True)
class FakeDownloader(Downloader):
    """
    Simulate a download that fails outright once,
    then returns success with the wrong file, then succeeds.
    """

    num_calls: int = 0

    def download(self, url: str, local_path: Path) -> bool:
        self.num_calls += 1
        if self.num_calls == 1:
            return False
        if self.num_calls == 2:
            local_path.write_text("AB")
            return True
        local_path.write_text(
            _TRUE_FILE_CONTENTS,
        )
        return True


def test_robust_downloader(tmp_path: Path) -> None:
    """
    test that the retry logic of the robust downloader works
    """
    dummy_file = tmp_path / "dummy"
    dummy_file.write_text(_TRUE_FILE_CONTENTS)
    out_path = tmp_path / "out"
    expected_hash = calc_md5_checksum(dummy_file)
    dl = FakeDownloader()
    robust_download(
        expected_hash,
        dest=out_path,
        url="fake_url",
        downloader=dl,
        max_backoff_time=0,
    )
    assert calc_md5_checksum(out_path) == expected_hash
    assert dl.num_calls == 3


@attrs.define(slots=True)
class AlwaysWrongFileDownloader(Downloader):
    """
    Simulate a server that reports success while returning the wrong bytes.
    """

    num_calls: int = 0

    def download(self, url: str, local_path: Path) -> bool:
        self.num_calls += 1
        local_path.write_text("wrong")
        return True


@attrs.define(slots=True)
class SkipsCompleteFileDownloader(Downloader):
    """
    Simulate aria2 with --continue=true: a file already at the full size is treated
    as finished, and is kept without being fetched again.

    The first call writes a corrupt file of the correct size.
    """

    num_calls: int = 0

    def download(self, url: str, local_path: Path) -> bool:
        self.num_calls += 1
        if local_path.exists() and local_path.stat().st_size == len(
            _TRUE_FILE_CONTENTS
        ):
            return True
        contents = "XYZ" if self.num_calls == 1 else _TRUE_FILE_CONTENTS
        local_path.write_text(contents)
        return True


@attrs.define(slots=True)
class NoFileThenGoodFileDownloader(Downloader):
    """
    Simulate a download that reports success without producing a file, then succeeds.
    """

    num_calls: int = 0

    def download(self, url: str, local_path: Path) -> bool:
        self.num_calls += 1
        if self.num_calls > 1:
            local_path.write_text(_TRUE_FILE_CONTENTS)
        return True


def test_robust_download_backs_off_after_hash_mismatch(tmp_path: Path) -> None:
    """
    Hash mismatches must wait between attempts, or a short outage exhausts every retry.
    """
    sleeps: list[float] = []
    dl = AlwaysWrongFileDownloader()
    with pytest.raises(DownloadAttemptFailedError):
        robust_download(
            md5sum="not_the_real_hash",
            dest=tmp_path / "out",
            url="fake_url",
            downloader=dl,
            max_outer_retries=3,
            sleep=sleeps.append,
        )
    assert dl.num_calls == 3
    assert len(sleeps) == 2


def test_robust_download_refetches_after_same_size_corrupt_file(
    tmp_path: Path,
) -> None:
    """
    A rejected file must not survive to the next attempt, where aria2 would keep it.
    """
    out_path = tmp_path / "out"
    robust_download(
        md5sum=_TRUE_FILE_MD5,
        dest=out_path,
        url="fake_url",
        downloader=SkipsCompleteFileDownloader(),
        max_outer_retries=3,
        sleep=_no_sleep,
    )
    assert out_path.read_text() == _TRUE_FILE_CONTENTS


def test_robust_download_retries_when_reported_success_left_no_file(
    tmp_path: Path,
) -> None:
    """
    Inspecting the rejected file must tolerate there being no file.
    """
    out_path = tmp_path / "out"
    robust_download(
        md5sum=_TRUE_FILE_MD5,
        dest=out_path,
        url="fake_url",
        downloader=NoFileThenGoodFileDownloader(),
        max_outer_retries=3,
        sleep=_no_sleep,
    )
    assert out_path.read_text() == _TRUE_FILE_CONTENTS


def _flag_value(cmd: list[str], flag: str) -> str:
    return cmd[cmd.index(flag) + 1]


def test_aria_downloader_sets_connections_and_split(tmp_path: Path) -> None:
    """
    aria2 opens min(-x, -s) connections to a server, so both must be set to the requested count.
    """
    commands: list[list[str]] = []

    def recording_executor(cmd: list[str]) -> str:
        commands.append(cmd)
        return ""

    num_connections = 8
    dl = AriaDownloader(num_simil=num_connections, executor=recording_executor)
    assert dl.download("fake_url", local_path=tmp_path / "out")
    (cmd,) = commands
    assert _flag_value(cmd, "-x") == str(num_connections)
    assert _flag_value(cmd, "-s") == str(num_connections)


def test_aria_downloader_rejects_too_many_connections() -> None:
    with pytest.raises(AssertionError):
        AriaDownloader(num_simil=ARIA2_MAX_CONNECTIONS_PER_SERVER + 1)
