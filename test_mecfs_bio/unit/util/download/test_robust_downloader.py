from pathlib import Path

import attrs

from mecfs_bio.util.download.robust_download import Downloader, robust_download
from mecfs_bio.util.download.verify import calc_md5_checksum

_TRUE_FILE_CONTENTS = "ABC"


@attrs.define
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


def test_robust_downloader(tmp_path: Path):
    """
    test that the retry logic of the robust downloader works
    """
    dummy_file = tmp_path / "dummy"
    dummy_file.write_text(_TRUE_FILE_CONTENTS)
    out_path = tmp_path / "out"
    expected_hash = calc_md5_checksum(dummy_file)
    robust_download(
        expected_hash,
        dest=out_path,
        url="fake_url",
        downloader=FakeDownloader(),
        max_backoff_time=0,
    )
