"""
Do deCODE SomaScan files share one variant row order? Range-read the first 24 MB (compressed)
of several files across both SomaScan groups, decompress the prefix, and compare the Name
column (chr:pos:ref:alt) row by row over the common prefix. Also compares object sizes.

Usage: AWS_* creds in env; pixi exec --spec python=3.12 --spec boto3 -- python -I \
    check_row_order_prefix.py <listing_smp.txt> <listing_raw.txt>
(the listings are `aws s3 ls --recursive` outputs of the two prefixes)
"""

import sys
import zlib

import boto3
from botocore.config import Config

B = "largescaleplasma-2023"
s3 = boto3.client("s3", endpoint_url="https://s3-ext.decode.is", region_name="us-east-1",
                  config=Config(signature_version="s3v4", s3={"addressing_style": "path"}))


def keys(listing: str, picks: list[int]) -> list[str]:
    rows = [ln.split()[3] for ln in open(listing) if ln.strip() and not ln.lstrip().startswith("Total")]
    return [rows[i] for i in picks]


def names_prefix(key: str, nbytes: int = 24 * 1024 * 1024) -> list[bytes]:
    body = s3.get_object(Bucket=B, Key=key, Range=f"bytes=0-{nbytes - 1}")["Body"].read()
    text = zlib.decompressobj(31).decompress(body)
    lines = text.split(b"\n")[1:-1]  # drop header and the partial last line
    return [ln.split(b"\t", 3)[2] for ln in lines]


def main() -> None:
    smp, raw = sys.argv[1], sys.argv[2]
    picked = keys(smp, [0, 40, 2000, 5283]) + keys(raw, [0, 40, 3000])
    ref_key, ref = picked[0], names_prefix(picked[0])
    print(f"reference {ref_key}: {len(ref):,} rows in prefix")
    for k in picked[1:]:
        other = names_prefix(k)
        n = min(len(ref), len(other))
        first_diff = next((i for i in range(n) if ref[i] != other[i]), None)
        print(f"{k.split('/')[-1][:60]:60s} rows={len(other):,} common={n:,} "
              f"first_mismatch={first_diff}")


if __name__ == "__main__":
    main()
