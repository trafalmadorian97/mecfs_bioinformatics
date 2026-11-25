"""
One sequential rapidgzip pass over a downloaded deCODE protein file, producing everything the
per-protein task stores:

  slim: Beta, SE and N of the shared index's rows, in index order, matched by Name (left join
      of the index's NAME on the file's Name). Protein files are not row-for-row copies of the
      index's template: some lack rows. An index row the file lacks gets nulls. A Name occurring
      twice among the matched rows fails the pass (the join is validated one-to-one), and so
      does a matched row whose effectAllele/otherAllele differ from the index's EA/NEA.
      To avoid parsing all ~33M lines, each block first parses only its Pos field and keeps the
      lines at a position of some index row (any chromosome); only those are parsed fully.
  posmap: the first line starting at or after every posmap_step bytes of uncompressed text, with
      its chromosome code, position, uncompressed offset and row number. Window readers use it
      to find the uncompressed span of a genomic window.
  significant_rows: every row, on any chromosome and in or out of the index, with
      minus_log10_pval >= DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P. The screen parses only the
      minus_log10_pval field (together with Pos, in one read); matching lines are parsed fully.
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
    DECODE_EFFECT_ALLELE_COL,
    DECODE_HEADER,
    DECODE_INDEX_NAME_COL,
    DECODE_MLOG10P_COL,
    DECODE_N_COL,
    DECODE_NAME_COL,
    DECODE_OTHER_ALLELE_COL,
    DECODE_POS_COL,
    DECODE_PVAL_COL,
    DECODE_RSIDS_COL,
    DECODE_SE_COL,
    DECODE_SIGNIFICANT_ROWS_MIN_MLOG10P,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
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
_MATCHED_COL = "__matched__"
SIGNIFICANT_ROWS_COLUMNS = [
    c for c in DECODE_HEADER if c not in (DECODE_RSIDS_COL, DECODE_PVAL_COL)
] + [SIGNIFICANT_ROW_NUMBER_COL]
_POS_FIELD = DECODE_HEADER.index(DECODE_POS_COL)
_MLOG10P_FIELD = DECODE_HEADER.index(DECODE_MLOG10P_COL)
_SLIM_SOURCE_COLUMNS = [DECODE_BETA_COL, DECODE_SE_COL, DECODE_N_COL]


@frozen(slots=True)
class OnePassResult:
    """
    slim: Beta, SE, N (deCODE names, Float64) of the index rows, in index order; null for an
        index row the file lacks.
    index_rows_missing: the number of index rows the file lacks.
    rows_scanned: the file's data row count.
    header_end_offset: uncompressed offset of the first data line.
    """

    slim: pl.DataFrame
    index_rows_missing: int
    rows_scanned: int
    posmap: pl.DataFrame
    significant_rows: pl.DataFrame
    header_end_offset: int
    seek_index: bytes


@frozen(slots=True)
class _BlockOutputs:
    candidates: pl.DataFrame
    posmap: list[tuple[int, int, int, int]]
    significant: pl.DataFrame
    next_mark: int


@frozen(slots=True)
class _BlockScreen:
    """Per line of a block: its Pos and minus_log10_pval (null for ".")."""

    pos: np.ndarray
    mlog10p: pl.Series


def one_pass(
    gz_path: Path,
    index: pl.DataFrame,
    parallelization: int = 8,
    chunk_size: int = 16 * MIB,
    posmap_step: int = 1 * MIB,
) -> OnePassResult:
    """index: the shared variant index's POS, EA, NEA and NAME (NAME unique)."""
    index_positions = np.unique(index[GWASLAB_POS_COL].cast(pl.Int64).to_numpy())
    outputs: list[_BlockOutputs] = []
    rows = 0
    with rapidgzip.RapidgzipFile(
        str(gz_path), parallelization=parallelization, chunk_size=chunk_size
    ) as f:
        stream = line_blocks(f)
        next_mark = stream.header_end
        for block in stream.blocks:
            out = _process_block(block, index_positions, next_mark, posmap_step)
            outputs.append(out)
            next_mark = out.next_mark
            rows = block.first_row + block.n_rows
        seek_index = export_index_from_open_file(f)
    slim = match_index_rows(index, pl.concat([o.candidates for o in outputs]))
    return OnePassResult(
        slim=slim.select(_SLIM_SOURCE_COLUMNS),
        index_rows_missing=int(slim[_MATCHED_COL].not_().sum()),
        rows_scanned=rows,
        posmap=pl.DataFrame(
            [r for o in outputs for r in o.posmap], schema=POSMAP_SCHEMA, orient="row"
        ),
        significant_rows=pl.concat([o.significant for o in outputs]),
        header_end_offset=stream.header_end,
        seek_index=seek_index,
    )


def match_index_rows(index: pl.DataFrame, candidates: pl.DataFrame) -> pl.DataFrame:
    """
    For each index NAME, in index order, the Beta, SE and N of the candidate row with that
    Name (nulls when there is none), and whether there was one. Fails when a NAME matches two
    candidate rows, or a row whose alleles are not the index's EA (effect) and NEA.
    """
    names = index[DECODE_INDEX_NAME_COL]
    matched = candidates.filter(pl.col(DECODE_NAME_COL).is_in(names.implode())).select(
        DECODE_NAME_COL,
        DECODE_EFFECT_ALLELE_COL,
        DECODE_OTHER_ALLELE_COL,
        *_SLIM_SOURCE_COLUMNS,
        pl.lit(True).alias(_MATCHED_COL),
    )
    joined = (
        index.select(
            DECODE_INDEX_NAME_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        )
        .join(
            matched,
            left_on=DECODE_INDEX_NAME_COL,
            right_on=DECODE_NAME_COL,
            how="left",
            validate="1:1",
            maintain_order="left",
        )
        .with_columns(pl.col(_MATCHED_COL).fill_null(False))
    )
    disagreeing = joined.filter(
        pl.col(_MATCHED_COL)
        & (
            (pl.col(DECODE_EFFECT_ALLELE_COL) != pl.col(GWASLAB_EFFECT_ALLELE_COL))
            | (pl.col(DECODE_OTHER_ALLELE_COL) != pl.col(GWASLAB_NON_EFFECT_ALLELE_COL))
        )
    ).height
    assert disagreeing == 0, (
        f"{disagreeing} rows match an index Name but not its effect and other alleles"
    )
    return joined


def _process_block(
    block: LineBlock, index_positions: np.ndarray, next_mark: int, posmap_step: int
) -> _BlockOutputs:
    screen = _screen_block(block)
    posmap, next_mark = posmap_rows(block, next_mark, posmap_step)
    at_index_position = _rows_at_positions(block, screen.pos, index_positions)
    return _BlockOutputs(
        candidates=parse_decode_lines(gather_lines(block, at_index_position)),
        posmap=posmap,
        significant=significant_rows(block, screen.mlog10p),
        next_mark=next_mark,
    )


def _screen_block(block: LineBlock) -> _BlockScreen:
    fields = pl.read_csv(
        io.BytesIO(block.text),
        separator="\t",
        has_header=False,
        columns=[_POS_FIELD, _MLOG10P_FIELD],
        infer_schema=False,
        null_values=["."],
        quote_char=None,
    )
    assert fields.height == block.n_rows
    return _BlockScreen(
        pos=fields.to_series(0).cast(pl.Int64).to_numpy(),
        mlog10p=fields.to_series(1).cast(pl.Float64),
    )


def _rows_at_positions(
    block: LineBlock, pos: np.ndarray, positions: np.ndarray
) -> np.ndarray:
    """Row numbers of the block's lines whose Pos is in positions (sorted, unique)."""
    if positions.size == 0:
        return np.empty(0, dtype=np.int64)
    at = np.minimum(np.searchsorted(positions, pos), positions.size - 1)
    return np.flatnonzero(positions[at] == pos) + block.first_row


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


def significant_rows(block: LineBlock, mlog10p: pl.Series) -> pl.DataFrame:
    """
    The block's rows with minus_log10_pval (mlog10p, one per line) at or above the threshold,
    in file order.
    """
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
