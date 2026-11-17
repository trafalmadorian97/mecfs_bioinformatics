"""
Interpret the output of an rsync copy.

rsync exits with code 23, partial transfer, for any per-file problem: an unreadable source,
an I/O error, or merely a failure to set permissions on the destination.  The last is
harmless when the destination is a mount that refuses chmod, such as a Windows drive
mounted in WSL without the metadata option.  rsync 3.5.x reports those chmod failures even
under --no-perms, where 3.4.x ignored them, so every copy to such a mount ends in code 23.
The functions here distinguish that case from real partial transfers, so callers can
tolerate the first without masking the second.
"""

import re

RSYNC_PARTIAL_TRANSFER_EXIT_CODE = 23

# rsync prefixes each per-file problem with "rsync: ", optionally followed by the role of
# the process that hit it, such as [receiver].  Progress output written with carriage
# returns can precede the message on the same line, so the match is not anchored to the
# start of a line.  The closing summary starts with "rsync error:" and is not matched.
_RSYNC_MESSAGE = re.compile(r"rsync: [^\n]*")
_PERMISSION_FAILURE = re.compile(
    r"rsync: (\[\w+\] )?failed to set permissions on .*: "
    r"Operation not permitted \(1\)"
)


def rsync_messages(output: str) -> list[str]:
    """
    Return the per-file problem messages in the combined output of an rsync run, one per
    message, with surrounding whitespace removed.
    """
    return [message.strip() for message in _RSYNC_MESSAGE.findall(output)]


def is_permission_only_partial_transfer(exit_code: int, output: str) -> bool:
    """
    Return True if an rsync run ended in a partial transfer solely because it could not
    set permissions on the destination.  Every file's data was still written, so a caller
    may treat the copy as complete, subject to its own verification.

    Any other exit code, or any problem message other than a permission failure, gives
    False, so a real partial transfer is never mistaken for this harmless one.
    """
    if exit_code != RSYNC_PARTIAL_TRANSFER_EXIT_CODE:
        return False
    messages = rsync_messages(output)
    return bool(messages) and all(
        _PERMISSION_FAILURE.fullmatch(message) for message in messages
    )
