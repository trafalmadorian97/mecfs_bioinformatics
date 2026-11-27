"""
Cut a sequential stream of deCODE sumstats text into blocks of whole data lines.

The stream is read in large pieces (READ_BYTES). Each yielded LineBlock holds whole lines only,
with its offset in the uncompressed text, the row number of its first data line (0 = the first
line after the header), and the position of every line's terminating newline, found with one
numpy scan. The header is checked and skipped.
"""

from collections.abc import Iterator
from typing import BinaryIO

import numpy as np
from attrs import field, frozen

from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    check_header,
)

MIB = 1024 * 1024
READ_BYTES = 16 * MIB
_NEWLINE = 10


@frozen(slots=True)
class LineBlock:
    """
    text: whole lines, each ending in a newline.
    offset: offset of text[0] in the uncompressed file.
    first_row: data row number of the first line.
    line_ends: positions in text of each line's newline (length = number of lines).
    """

    text: bytes
    offset: int
    first_row: int
    line_ends: np.ndarray = field(eq=False)

    def __attrs_post_init__(self) -> None:
        assert self.line_ends.ndim == 1 and self.line_ends.dtype.kind == "i"

    @property
    def n_rows(self) -> int:
        return len(self.line_ends)

    @property
    def end_offset(self) -> int:
        return self.offset + len(self.text)

    @property
    def line_starts(self) -> np.ndarray:
        return np.concatenate(([0], self.line_ends[:-1] + 1))


@frozen(slots=True)
class HeaderedStream:
    """header_end: offset of the first data line (just past the header's newline)."""

    header_end: int
    blocks: Iterator[LineBlock]


def _blocks(
    f: BinaryIO, pending: bytes, offset: int, read_bytes: int
) -> Iterator[LineBlock]:
    row = 0
    eof = False
    while not eof:
        more = f.read(read_bytes)
        eof = not more
        pending += more
        if not eof and len(pending) < read_bytes:
            continue
        if eof and pending and not pending.endswith(b"\n"):
            pending += b"\n"
        cut = len(pending) if eof else pending.rindex(b"\n") + 1
        text, pending = pending[:cut], pending[cut:]
        if text:
            ends = np.flatnonzero(np.frombuffer(text, dtype=np.uint8) == _NEWLINE)
            yield LineBlock(text=text, offset=offset, first_row=row, line_ends=ends)
            offset += len(text)
            row += len(ends)


def line_blocks(f: BinaryIO, read_bytes: int = READ_BYTES) -> HeaderedStream:
    """Check the header of the decompressed stream f, then yield its data lines in blocks."""
    first = f.read(read_bytes)
    header_end = first.index(b"\n") + 1
    check_header(first[: header_end - 1])
    return HeaderedStream(
        header_end=header_end,
        blocks=_blocks(f, first[header_end:], header_end, read_bytes),
    )


def gather_lines(block: LineBlock, rows: np.ndarray) -> bytes:
    """The text of the given data row numbers (ascending, all inside block), concatenated."""
    local = rows.astype(np.int64) - block.first_row
    ends = block.line_ends[local]
    starts = np.where(local == 0, 0, block.line_ends[np.maximum(local - 1, 0)] + 1)
    return b"".join(
        block.text[s : e + 1]
        for s, e in zip(starts.tolist(), ends.tolist(), strict=True)
    )
