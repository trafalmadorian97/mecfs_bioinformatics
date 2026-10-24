from pathlib import Path

import attrs
import pytest

from mecfs_bio.util.download.robust_download import (
    ARIA2_MAX_CONNECTIONS_PER_SERVER,
    AriaDownloader,
    Downloader,
    robust_download,
)
from mecfs_bio.util.download.verify import calc_md5_checksum

_TRUE_FILE_CONTENTS = "ABC"


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
