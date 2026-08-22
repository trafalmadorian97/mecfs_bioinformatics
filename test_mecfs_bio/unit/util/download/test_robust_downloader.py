from pathlib import Path

import attrs
import pytest

from mecfs_bio.util.download.robust_download import (
    Downloader,
    aria_command,
    robust_download,
)
from mecfs_bio.util.download.verify import calc_md5_checksum

_TRUE_FILE_CONTENTS = "ABC"


@attrs.define
class FakeDownloader(Downloader):
    """
    Simulate a download that fails outright once,
    then returns success with the wrong file, then succeeds.
    """

    num_calls: int = 0

    def download(
        self, url: str, local_path: Path, request_connections: int | None = None
    ) -> bool:
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


def test_robust_downloader(tmp_path: Path):
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


def test_aria_command_connections_set_both_x_and_s(tmp_path: Path):
    # aria2 caps the real connection count at min(-x, -s), so the requested count must
    # appear on both flags for a multi-connection download to actually happen.
    cmd = aria_command(
        url="http://host/f",
        local_path=tmp_path / "f",
        connections=16,
        summary_interval=10,
    )
    assert cmd[cmd.index("-x") + 1] == "16"
    assert cmd[cmd.index("-s") + 1] == "16"


def test_aria_command_rejects_too_many_connections(tmp_path: Path):
    # aria2 caps connections per server at 16; a larger request must fail fast here
    # rather than erroring mid-download.
    with pytest.raises(AssertionError):
        aria_command(
            url="http://host/f",
            local_path=tmp_path / "f",
            connections=17,
            summary_interval=10,
        )
