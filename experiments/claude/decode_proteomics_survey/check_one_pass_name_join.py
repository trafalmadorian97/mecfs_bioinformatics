"""
The Name-join one_pass on the local AP4M1 reference file: time, peak RSS, index rows missing, and
slim == a brute-force left join of the index NAME on the file's Name.

  pixi r python experiments/claude/decode_proteomics_survey/check_one_pass_name_join.py
"""

import resource
import time

import polars as pl
from polars.testing import assert_frame_equal

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_index import (
    HAPMAP_3_DECODE_SMP_INDEX,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_reference_protein_file import (
    DECODE_SMP_AP4M1_RAW_FILE,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.decode_database.line_blocks import MIB
from mecfs_bio.build_system.task.decode_database.one_pass import one_pass
from mecfs_bio.constants.decode_somascan_constants import DECODE_INDEX_NAME_COL
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL


def main() -> None:
    assets = DEFAULT_RUNNER.run([DECODE_SMP_AP4M1_RAW_FILE, HAPMAP_3_DECODE_SMP_INDEX])
    raw = assets[DECODE_SMP_AP4M1_RAW_FILE.asset_id]
    index_asset = assets[HAPMAP_3_DECODE_SMP_INDEX.asset_id]
    assert isinstance(raw, FileAsset) and isinstance(index_asset, FileAsset)
    index = pl.read_parquet(index_asset.path, columns=[GWASLAB_POS_COL, DECODE_INDEX_NAME_COL])
    t0 = time.perf_counter()
    result = one_pass(raw.path, index, chunk_size=4 * MIB)
    elapsed = time.perf_counter() - t0
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(
        f"one_pass {elapsed:.1f}s, peak RSS {peak:.0f} MB, rows {result.rows_scanned}, "
        f"index rows missing {result.index_rows_missing}, significant {result.significant_rows.height}",
        flush=True,
    )
    full = pl.read_csv(
        raw.path,
        separator="\t",
        columns=["Name", "Beta", "SE", "N"],
        null_values=["."],
        quote_char=None,
    )
    expected = index.join(
        full.filter(pl.col("Name").is_in(index[DECODE_INDEX_NAME_COL].implode())),
        left_on=DECODE_INDEX_NAME_COL,
        right_on="Name",
        how="left",
        maintain_order="left",
    ).select("Beta", "SE", pl.col("N").cast(pl.Float64))
    assert_frame_equal(result.slim, expected)
    print(f"slim == brute force ({result.slim.height} rows)")


if __name__ == "__main__":
    main()
