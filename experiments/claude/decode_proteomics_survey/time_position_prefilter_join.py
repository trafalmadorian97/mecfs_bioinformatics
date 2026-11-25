"""
Matching index rows by Name with a cheap position pre-filter, inside the per-protein pass:
per block, parse only Pos, keep rows at any index position (on any chromosome), gather and fully parse
those lines (as the current gather does), then left-join the index on Name with validate 1:1.
Timed on the local AP4M1 copy against the decompress-and-block pass alone.

  pixi r python experiments/claude/decode_proteomics_survey/time_position_prefilter_join.py
"""

import io
import resource
import time

import numpy as np
import polars as pl
import rapidgzip

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.decode_sumstats_format import (
    parse_decode_lines,
)
from mecfs_bio.build_system.task.decode_database.line_blocks import gather_lines, line_blocks



def peak_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX])
    raw, index_asset = (
        assets[t.asset_id] for t in (DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX)
    )
    assert isinstance(raw, FileAsset) and isinstance(index_asset, FileAsset)
    index = pl.read_parquet(index_asset.path, columns=["CHR", "POS", "NAME", "SOURCE_ROW"])
    keys = np.unique(index["POS"].cast(pl.Int64).to_numpy())

    t0 = time.perf_counter()
    screen = gather = 0.0
    kept = []
    with rapidgzip.RapidgzipFile(str(raw.path), parallelization=8) as f:
        for block in line_blocks(f).blocks:
            s0 = time.perf_counter()
            row_keys = (
                pl.read_csv(
                    io.BytesIO(block.text),
                    separator="\t",
                    has_header=False,
                    columns=[1],
                    infer_schema=False,
                    quote_char=None,
                )
                .to_series(0)
                .cast(pl.Int64)
                .to_numpy()
            )
            at = np.minimum(np.searchsorted(keys, row_keys), len(keys) - 1)
            rows = np.flatnonzero(keys[at] == row_keys) + block.first_row
            s1 = time.perf_counter()
            parsed = parse_decode_lines(gather_lines(block, rows)).with_columns(
                pl.Series("file_row", rows, dtype=pl.UInt64)
            )
            kept.append(parsed.select("Name", "Beta", "SE", "N", "file_row"))
            screen += s1 - s0
            gather += time.perf_counter() - s1
    j0 = time.perf_counter()
    candidates = pl.concat(kept)
    joined = index.join(
        candidates.filter(pl.col("Name").is_in(index["NAME"].implode())),
        left_on="NAME",
        right_on="Name",
        how="left",
        validate="1:1",
    )
    join = time.perf_counter() - j0
    total = time.perf_counter() - t0
    print(
        f"pass with position pre-filter + join: {total:.1f}s; Chrom/Pos screen {screen:.1f}s, "
        f"gather+parse of {candidates.height} candidate rows {gather:.1f}s, join {join:.2f}s "
        f"(peak RSS {peak_mb():.0f} MB)"
    )
    print(
        f"index rows {index.height}, matched {joined['file_row'].is_not_null().sum()}, "
        f"in file order: {joined['file_row'].drop_nulls().is_sorted()}, "
        f"file_row == SOURCE_ROW for all: {(joined['file_row'] == joined['SOURCE_ROW']).all()}"
    )


if __name__ == "__main__":
    main()
