"""
Measure how much slower rsID assignment gets after dbSNP moved to the external drive.

The rsID assignment join (annovar_37_basic_rsid_assignment) reads five columns of the
353M-row annovar dbSNP150 parquet through a lazy duckdb scan.  After the asset store
split, that file lives on /mnt/d instead of local ext4, so the question is how much wall
clock the extra I/O costs.

Raw dd overstates the answer, because parquet decode is CPU work that overlaps the read.
This script measures the real scan twice:

    cold pass -- data comes off /mnt/d, so the scan is bounded by whichever of I/O and
                 decode is slower
    warm pass -- the same bytes are served from page cache, which is the closest available
                 proxy for the old local-disk behaviour now that the local copy is gone

The difference between the two passes is the honest estimate of the regression.
"""

import subprocess
import time
from pathlib import Path

import duckdb

DBSNP_PARQUET = Path(
    "/mnt/d/asset_store_remote/reference_data/db_snp_reference_data/build_37/annovar/"
    "db_snp150_annovar_proc_parquet_unique_direct_download.parquet"
)

# The columns the join actually touches: right_on plus the rsid it carries through.
JOIN_COLUMNS = ("int_chrom", "POS", "ALT", "REF", "rsid")


def scan_join_columns(path: Path) -> tuple[int, float]:
    """
    Read exactly the columns the join reads, forcing full materialization of each, and
    return the row count together with the elapsed seconds.

    The aggregate must depend on the column VALUES.  count(column) is answered from the
    parquet footer statistics without reading a single data page, which makes the query
    finish in a fraction of a second and measures nothing.  Hashing every join key forces
    the scan the real join performs.
    """
    keys = ", ".join(JOIN_COLUMNS)
    query = (
        f"select count(*), sum(hash({keys})::HUGEINT) from read_parquet('{path}')"
    )
    start = time.monotonic()
    result = duckdb.sql(query).fetchone()
    elapsed = time.monotonic() - start
    assert result is not None
    return result[0], elapsed


def drop_page_cache_for(path: Path) -> None:
    """
    Evict the file from the page cache so the next read is genuinely cold.  'dd count=0
    iflag=nocache' applies POSIX_FADV_DONTNEED to the whole file rather than to a single
    block.
    """
    subprocess.run(["sync"], capture_output=True, check=False)
    subprocess.run(
        ["dd", f"if={path}", "of=/dev/null", "count=0", "iflag=nocache"],
        capture_output=True,
        check=False,
    )


def main() -> None:
    size_gb = DBSNP_PARQUET.stat().st_size / 1e9
    print(f"file: {DBSNP_PARQUET}")
    print(f"size: {size_gb:.2f} GB")
    print(f"columns scanned: {', '.join(JOIN_COLUMNS)}\n")

    drop_page_cache_for(DBSNP_PARQUET)
    rows, cold = scan_join_columns(DBSNP_PARQUET)
    print(f"cold pass (served from /mnt/d):     {cold:7.1f} s   rows={rows:,}")

    _, warm = scan_join_columns(DBSNP_PARQUET)
    print(f"warm pass (served from page cache): {warm:7.1f} s")

    print(
        f"\nestimated regression per rsID assignment run: "
        f"+{cold - warm:.1f} s ({(cold - warm) / 60:.1f} min)"
    )
    print(f"cold/warm ratio: {cold / warm:.1f}x")


if __name__ == "__main__":
    main()
