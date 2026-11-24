"""
Which deCODE SomaScan files differ in uncompressed size from the group's template? Reads each
file's gzip trailer (ISIZE: uncompressed length mod 2^32, last 4 bytes) with one range request.

  pixi r python experiments/claude/decode_proteomics_survey/isize_survey.py <smp|raw> <out.csv>
"""

import sys
from concurrent.futures import ThreadPoolExecutor

import polars as pl

from mecfs_bio.asset_generator.decode_somascan_manifest import load_decode_somascan_proteins
from mecfs_bio.build_system.wf.s3_access import Boto3S3Access
from mecfs_bio.util.download.parallel_ranges import (
    REQUESTS_RANGE_GETTER,
    fetch_range_with_retries,
)


def main() -> None:
    group, out = sys.argv[1], sys.argv[2]
    proteins = load_decode_somascan_proteins(group)  # ty: ignore[invalid-argument-type]
    access = Boto3S3Access()

    def isize(p) -> int:  # noqa: ANN001
        url = access.presigned_url(p.ref, 3600)
        tail = fetch_range_with_retries(REQUESTS_RANGE_GETTER, url, p.size - 4, p.size - 1)
        return int.from_bytes(tail, "little")

    with ThreadPoolExecutor(32) as pool:
        sizes = list(pool.map(isize, proteins))
    table = pl.DataFrame(
        {
            "seq_id": [p.seq_id for p in proteins],
            "gene": [p.gene_symbol for p in proteins],
            "isize": sizes,
        }
    )
    table.write_csv(out)
    template = sizes[0]
    differ = table.filter(pl.col("isize") != template)
    print(f"{group}: {table.height} files, template isize {template}, {differ.height} differ")
    print(differ.head(20))


if __name__ == "__main__":
    main()
