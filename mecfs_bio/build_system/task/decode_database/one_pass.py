"""
One sequential rapidgzip pass over a downloaded deCODE protein file, producing everything the
per-protein task stores:

  slim: Beta, SE and N of the shared index's rows (gathered by row number, not by a join), in
      index order. Every gathered Name must equal the index's NAME and none may carry the
      missing-allele marker, and the file must have the index's total row count, so a file whose
      rows are reordered or differ in number fails instead of yielding misaligned rows.
  posmap: the first line starting at or after every posmap_step bytes of uncompressed text, with
      its chromosome code, position, uncompressed offset and row number. Window readers use it
      to find the uncompressed span of a genomic window.
  significant_rows: every row, on any chromosome and in or out of the index, with
      minus_log10_pval >= DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P. Only the minus_log10_pval
      column of each block is parsed for the screen; the matching lines are then parsed fully.
  seek_index: the zstd-compressed GZIDX index rapidgzip builds during the pass (one seek point
      per chunk_size of compressed input), importable by indexed_gzip.
"""

import io
from pathlib import Path

import numpy as np
import polars as pl
import rapidgzip
from attrs import frozen

from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    chrom_code,
    parse_decode_lines,
)
from mecfs_bio.build_system.task.decode_database.line_blocks import (
    MIB,
    LineBlock,
    gather_lines,
    line_blocks,
)
from mecfs_bio.constants.decode_somascan_constants import (
    DECODE_BETA_COL,
    DECODE_HEADER,
    DECODE_INDEX_NAME_COL,
    DECODE_INDEX_SOURCE_ROW_COL,
    DECODE_MISSING_ALLELE,
    DECODE_MLOG10P_COL,
    DECODE_N_COL,
    DECODE_NAME_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_PVAL_COL,
    DECODE_RSIDS_COL,
    DECODE_SE_COL,
    DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P,
)
from mecfs_bio.util.file_io.gzip_seek_index import export_index_from_open_file

POSMAP_CHROM_COL = "chrom"
POSMAP_POS_COL = "pos"
POSMAP_LINE_OFFSET_COL = "line_offset"
POSMAP_ROW_NUMBER_COL = "row_number"
POSMAP_SCHEMA = pl.Schema(
    {
        POSMAP_CHROM_COL: pl.Int8,
        POSMAP_POS_COL: pl.Int64,
        POSMAP_LINE_OFFSET_COL: pl.UInt64,
        POSMAP_ROW_NUMBER_COL: pl.UInt64,
    }
)
SIGNIFICANT_ROW_NUMBER_COL = "row_number"
SIGNIFICANT_ROWS_COLUMNS = [
    c for c in DECODE_HEADER if c not in (DECODE_RSIDS_COL, DECODE_PVAL_COL)
] + [SIGNIFICANT_ROW_NUMBER_COL]
_MLOG10P_FIELD = DECODE_HEADER.index(DECODE_MLOG10P_COL)
_SLIM_SOURCE_COLUMNS = [DECODE_BETA_COL, DECODE_SE_COL, DECODE_N_COL]


@frozen(slots=True)
class OnePassResult:
    """
    slim: Beta, SE, N (deCODE names, Float64) of the index rows, in index order.
    header_end_offset: uncompressed offset of the first data line.
    """

    slim: pl.DataFrame
    posmap: pl.DataFrame
    significant_rows: pl.DataFrame
    header_end_offset: int
    seek_index: bytes


@frozen(slots=True)
class _BlockOutputs:
    gathered: pl.DataFrame
    posmap: list[tuple[int, int, int, int]]
    significant: pl.DataFrame
    next_mark: int


def one_pass(
    gz_path: Path,
    index: pl.DataFrame,
    total_rows: int,
    parallelization: int = 8,
    chunk_size: int = 16 * MIB,
    posmap_step: int = 1 * MIB,
) -> OnePassResult:
    """index: the shared variant index (SOURCE_ROW ascending, NAME)."""
    targets = index[DECODE_INDEX_SOURCE_ROW_COL].to_numpy().astype(np.int64)
    assert np.all(np.diff(targets) > 0), "index SOURCE_ROW must be strictly ascending"
    outputs: list[_BlockOutputs] = []
    rows = 0
    with rapidgzip.RapidgzipFile(
        str(gz_path), parallelization=parallelization, chunk_size=chunk_size
    ) as f:
        stream = line_blocks(f)
        next_mark = stream.header_end
        for block in stream.blocks:
            out = _process_block(block, targets, next_mark, posmap_step)
            outputs.append(out)
            next_mark = out.next_mark
            rows = block.first_row + block.n_rows
        seek_index = export_index_from_open_file(f)
    assert rows == total_rows, (
        f"{gz_path.name} has {rows} rows, the index's template {total_rows}: row order differs"
    )
    gathered = pl.concat([o.gathered for o in outputs])
    _check_gathered(gathered, index)
    return OnePassResult(
        slim=gathered.select(_SLIM_SOURCE_COLUMNS),
        posmap=pl.DataFrame(
            [r for o in outputs for r in o.posmap], schema=POSMAP_SCHEMA, orient="row"
        ),
        significant_rows=pl.concat([o.significant for o in outputs]),
        header_end_offset=stream.header_end,
        seek_index=seek_index,
    )


def _process_block(
    block: LineBlock, targets: np.ndarray, next_mark: int, posmap_step: int
) -> _BlockOutputs:
    lo = np.searchsorted(targets, block.first_row)
    hi = np.searchsorted(targets, block.first_row + block.n_rows)
    posmap, next_mark = posmap_rows(block, next_mark, posmap_step)
    return _BlockOutputs(
        gathered=parse_decode_lines(gather_lines(block, targets[lo:hi])),
        posmap=posmap,
        significant=significant_rows(block),
        next_mark=next_mark,
    )


def posmap_rows(
    block: LineBlock, next_mark: int, posmap_step: int
) -> tuple[list[tuple[int, int, int, int]], int]:
    """
    Position-map rows (chrom code, pos, line offset, row number) for the first line starting at
    or after each mark inside block; returns them with the next mark.
    """
    starts = block.line_starts
    rows: list[tuple[int, int, int, int]] = []
    while next_mark < block.end_offset:
        i = int(np.searchsorted(starts, next_mark - block.offset, side="left"))
        if i >= len(starts):
            break
        s = int(starts[i])
        fields = block.text[s : s + 64].split(b"\t", 2)
        rows.append(
            (
                chrom_code(fields[0].decode()),
                int(fields[1]),
                block.offset + s,
                block.first_row + i,
            )
        )
        next_mark = block.offset + s + posmap_step
    return rows, next_mark


def significant_rows(block: LineBlock) -> pl.DataFrame:
    """The block's rows with minus_log10_pval at or above the threshold, in file order."""
    mlog10p = (
        pl.read_csv(
            io.BytesIO(block.text),
            separator="\t",
            has_header=False,
            columns=[_MLOG10P_FIELD],
            infer_schema=False,
            null_values=["."],
            quote_char=None,
        )
        .to_series(0)
        .cast(pl.Float64)
    )
    assert mlog10p.len() == block.n_rows
    local = np.flatnonzero(
        (mlog10p >= DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P).fill_null(False).to_numpy()
    )
    row_numbers = local + block.first_row
    return (
        parse_decode_lines(gather_lines(block, row_numbers))
        .with_columns(
            pl.Series(SIGNIFICANT_ROW_NUMBER_COL, row_numbers, dtype=pl.UInt64)
        )
        .select(SIGNIFICANT_ROWS_COLUMNS)
    )


def _check_gathered(gathered: pl.DataFrame, index: pl.DataFrame) -> None:
    assert gathered.height == index.height, (
        f"gathered {gathered.height} rows for an index of {index.height}"
    )
    mismatched = int((gathered[DECODE_NAME_COL] != index[DECODE_INDEX_NAME_COL]).sum())
    assert mismatched == 0, (
        f"{mismatched} gathered rows do not match the index's Name: row order differs"
    )
    missing = int((gathered[DECODE_OTHER_ALLELE_COL] == DECODE_MISSING_ALLELE).sum())
    assert missing == 0, f"{missing} gathered rows carry the missing-allele marker"
