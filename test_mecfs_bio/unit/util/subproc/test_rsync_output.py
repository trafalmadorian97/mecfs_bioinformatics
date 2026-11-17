import pytest

from mecfs_bio.util.subproc.rsync_output import (
    RSYNC_PARTIAL_TRANSFER_EXIT_CODE,
    is_permission_only_partial_transfer,
)

# Shape taken from a real rsync 3.5.2 copy to a WSL DrvFs mount: progress written with
# carriage returns runs straight into the message on the same line.
_PERMISSION_ONLY_OUTPUT = (
    'rsync: [generator] failed to set permissions on "/mnt/d/store/x/.": '
    "Operation not permitted (1)\n"
    "   952MB   9%   71.7MB/s    0:00:13 (xfr#1, to-chk=7/9)"
    'rsync: [receiver] failed to set permissions on "/mnt/d/store/x/.a.pickle.3rsYpU": '
    "Operation not permitted (1)\r\n"
    "   9.86GB 100%   67.3MB/s    0:02:26 (xfr#8, to-chk=0/9)\n"
    "rsync error: some files/attrs were not transferred (see previous errors) "
    "(code 23) at main.c(1416) [sender=3.5.2]\n"
)

_UNREADABLE_SOURCE = (
    'rsync: [sender] send_files failed to open "/src/x/b.pickle": '
    "Permission denied (13)\n"
)


def test_permission_failures_alone_are_tolerated() -> None:
    assert is_permission_only_partial_transfer(
        exit_code=RSYNC_PARTIAL_TRANSFER_EXIT_CODE, output=_PERMISSION_ONLY_OUTPUT
    )


@pytest.mark.parametrize(
    argnames=("exit_code", "output"),
    argvalues=[
        (
            RSYNC_PARTIAL_TRANSFER_EXIT_CODE,
            _PERMISSION_ONLY_OUTPUT + _UNREADABLE_SOURCE,
        ),
        (RSYNC_PARTIAL_TRANSFER_EXIT_CODE, ""),
        (11, _PERMISSION_ONLY_OUTPUT),
    ],
    ids=["another_error_alongside", "no_messages", "different_exit_code"],
)
def test_anything_else_is_a_real_failure(exit_code: int, output: str) -> None:
    assert not is_permission_only_partial_transfer(exit_code=exit_code, output=output)
