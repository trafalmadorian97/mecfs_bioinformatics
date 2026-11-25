"""
What would matching index rows by Name (a polars join) add to the per-protein pass? Times, on the
local AP4M1 copy, the extra per-block parsing of Name, Beta, SE and N for every row, and the
join with the HapMap3 index (left join from the index, after keeping index Names, validate 1:1),
against the decompress-and-block pass alone.

  pixi r python experiments/claude/decode_proteomics_survey/time_name_join.py
"""

import io
import resource
import time

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
from mecfs_bio.build_system.task.decode_database.line_blocks import line_blocks
from mecfs_bio.constants.decode_somascan_constants import DECODE_HEADER

FIELDS = ["Name", "Beta", "SE", "N"]
FIELD_INDICES = [DECODE_HEADER.index(f) for f in FIELDS]


def peak_mb() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX])
    raw, index_asset = (assets[t.asset_id] for t in (DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX))
    assert isinstance(raw, FileAsset) and isinstance(index_asset, FileAsset)
    index = pl.read_parquet(index_asset.path, columns=["NAME", "SOURCE_ROW"])
    names = index["NAME"].implode()

    t0 = time.perf_counter()
    with rapidgzip.RapidgzipFile(str(raw.path), parallelization=8) as f:
        for _ in line_blocks(f).blocks:
            pass
    baseline = time.perf_counter() - t0
    print(f"decompress-and-block pass alone: {baseline:.1f}s (peak RSS {peak_mb():.0f} MB)")

    t0 = time.perf_counter()
    parse = 0.0
    kept = []
    with rapidgzip.RapidgzipFile(str(raw.path), parallelization=8) as f:
        for block in line_blocks(f).blocks:
            p0 = time.perf_counter()
            frame = pl.read_csv(
                io.BytesIO(block.text),
                separator="\t",
                has_header=False,
                columns=FIELD_INDICES,
                infer_schema=False,
                quote_char=None,
            )
            frame = (
                frame.rename(dict(zip(frame.columns, FIELDS, strict=True)))
                .with_columns(
                    pl.col("Beta", "SE", "N").replace(".", None).cast(pl.Float64)
                )
                .with_row_index("file_row", offset=block.first_row)
            )
            kept.append(frame.filter(pl.col("Name").is_in(names)))
            parse += time.perf_counter() - p0
    j0 = time.perf_counter()
    matched = pl.concat(kept)
    joined = index.join(
        matched, left_on="NAME", right_on="Name", how="left", validate="1:1"
    )
    join = time.perf_counter() - j0
    total = time.perf_counter() - t0
    print(
        f"with Name parse + join: {total:.1f}s; per-block parse+filter {parse:.1f}s, "
        f"join {join:.2f}s (peak RSS {peak_mb():.0f} MB)"
    )
    print(
        f"index rows {index.height}, matched {joined['file_row'].is_not_null().sum()}, "
        f"in file order: {joined['file_row'].drop_nulls().is_sorted()}, "
        f"file_row == SOURCE_ROW for all: {(joined['file_row'] == joined['SOURCE_ROW']).all()}"
    )


if __name__ == "__main__":
    main()
