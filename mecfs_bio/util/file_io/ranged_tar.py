"""
Read individual members of an uncompressed tar archive through byte-range reads.

A tar archive is a sequence of members, each a 512-byte header followed by the member's data
padded up to a multiple of 512 bytes. The header records the member's size, so the position of
the next header is known without reading the data in between. That makes it possible to list
the members of a large remote tar with one small range read per member, and then fetch only
the members that are needed.

GNU tar stores names longer than 100 characters in an extra member of type L (named
././@LongLink) whose data is the full name of the member that follows it. Those are folded
into the following member. POSIX pax extended headers are not supported and fail fast.
"""

import tarfile
from typing import Protocol

from attrs import define, frozen

_TAR_BLOCK = tarfile.BLOCKSIZE
# Bytes fetched per hop between headers: one header plus a GNU long-name member and the
# header that follows it fit comfortably, so a long name costs no extra read.
_HOP_BYTES = 4096
_PAX_TYPES = (tarfile.XHDTYPE, tarfile.XGLTYPE, tarfile.SOLARIS_XHDTYPE)


class RangeReader(Protocol):
    def read_range(self, start: int, length: int) -> bytes:
        """Return the length bytes starting at byte offset start; fewer if the file ends first."""
        ...


@frozen(slots=True)
class TarMemberExtent:
    """Where a regular-file member's data lies: size bytes starting at data_offset."""

    name: str
    data_offset: int
    size: int


@define(slots=True)
class _HopReader:
    """
    Serves small reads near the current position from one cached range read.

    buffer holds bytes [buffer_start, buffer_start + len(buffer)) of the archive.
    """

    reader: RangeReader
    buffer_start: int = 0
    buffer: bytes = b""

    def read(self, offset: int, length: int) -> bytes:
        buffer_end = self.buffer_start + len(self.buffer)
        if not (self.buffer_start <= offset and offset + length <= buffer_end):
            self.buffer_start = offset
            self.buffer = self.reader.read_range(offset, max(length, _HOP_BYTES))
        chunk = self.buffer[
            offset - self.buffer_start : offset - self.buffer_start + length
        ]
        assert len(chunk) == length, (
            f"tar archive ended inside a member: wanted {length} bytes at offset {offset}"
        )
        return chunk


def _padded(size: int) -> int:
    return -(-size // _TAR_BLOCK) * _TAR_BLOCK


def list_tar_members(reader: RangeReader) -> tuple[TarMemberExtent, ...]:
    """
    The regular-file members of the tar archive behind reader, in archive order.

    Directories and other non-file members are skipped. Member data is never read, except
    the data of GNU long-name members, which holds the next member's name.
    """
    members: list[TarMemberExtent] = []
    hop_reader = _HopReader(reader=reader)
    offset = 0
    long_name: str | None = None
    while True:
        header = hop_reader.read(offset, _TAR_BLOCK)
        if header == b"\0" * _TAR_BLOCK:
            return tuple(members)
        info = tarfile.TarInfo.frombuf(header, "utf-8", "surrogateescape")
        assert info.type not in _PAX_TYPES, (
            f"pax extended headers are not supported (member at offset {offset})"
        )
        data_offset = offset + _TAR_BLOCK
        if info.type == tarfile.GNUTYPE_LONGNAME:
            name_bytes = hop_reader.read(data_offset, info.size)
            long_name = name_bytes.rstrip(b"\0").decode("utf-8", "surrogateescape")
        else:
            if info.isreg():
                members.append(
                    TarMemberExtent(
                        name=long_name if long_name is not None else info.name,
                        data_offset=data_offset,
                        size=info.size,
                    )
                )
            long_name = None
        offset = data_offset + _padded(info.size)


def read_tar_member(reader: RangeReader, member: TarMemberExtent) -> bytes:
    """The full data of one member, fetched with a single range read."""
    data = reader.read_range(member.data_offset, member.size)
    assert len(data) == member.size, (
        f"short read of tar member {member.name}: {len(data)} of {member.size} bytes"
    )
    return data
