"""
Probe whether a UKB-PPP protein tar on Synapse can be read one chromosome at a time.

Walks the tar headers with HTTP Range requests against the presigned S3 URL (one 512-byte
read per member, skipping member data), reports each member's name, data offset and size,
checks whether the gzip members are BGZF (blocked, hence seekable), and then range-reads one
chromosome member and decompresses it to confirm it is a complete regenie file.

Run: pixi r python experiments/claude/protein_coloc/probe_ppp_tar_range_layout.py 2>&1 | tee experiments/claude/protein_coloc/probe_ppp_tar_range_layout.log
"""

import gzip
import io
import tarfile
import time

import polars as pl
import requests
import synapseclient

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)

_BLOCK = 512
_PROBE_CHROM = "chr19_"


def _range(url: str, start: int, length: int) -> bytes:
    response = requests.get(
        url, headers={"Range": f"bytes={start}-{start + length - 1}"}, timeout=120
    )
    response.raise_for_status()
    assert response.status_code == 206, response.status_code
    return response.content


def _is_bgzf(gzip_head: bytes) -> bool:
    # BGZF: gzip with FLG.FEXTRA set and an extra subfield with id "BC".
    has_extra = gzip_head[3] & 0x04 != 0
    return has_extra and gzip_head[12:14] == b"BC"


def main() -> None:
    protein = load_ppp_protein_files(EUR_DISCOVERY_PPP_MANIFEST_PATH)[0]
    print(f"protein {protein.oid} synid {protein.synid}")
    client = synapseclient.login()
    entity = client.get(protein.synid, downloadFile=False)
    info = client._getFileHandleDownload(entity.dataFileHandleId, protein.synid)
    url = info["preSignedURL"]
    print(f"tar size {entity.properties.get('contentSize') or 'unknown'}")

    start_time = time.monotonic()
    offset = 0
    n_requests = 0
    members: list[tuple[str, int, int]] = []
    while True:
        header = _range(url, offset, _BLOCK)
        n_requests += 1
        if header == b"\0" * _BLOCK:
            break
        info_member = tarfile.TarInfo.frombuf(header, "utf-8", "surrogateescape")
        data_offset = offset + _BLOCK
        members.append((info_member.name, data_offset, info_member.size))
        offset = data_offset + -(-info_member.size // _BLOCK) * _BLOCK
    print(f"walked {len(members)} members in {n_requests} requests, "
          f"{time.monotonic() - start_time:.1f}s")
    total = sum(size for _, _, size in members)
    for name, data_offset, size in members:
        print(f"  {data_offset:>12} {size:>12} {name}")
    print(f"total member bytes {total}")

    regular = [m for m in members if m[2] > 0]
    gzip_head = _range(url, regular[0][1], 32)
    print(f"first member BGZF: {_is_bgzf(gzip_head)}  head {gzip_head[:16].hex()}")

    name, data_offset, size = next(m for m in regular if _PROBE_CHROM in m[0])
    start_time = time.monotonic()
    payload = _range(url, data_offset, size)
    print(f"fetched {name} ({size} bytes) in {time.monotonic() - start_time:.1f}s")
    frame = pl.read_csv(io.BytesIO(gzip.decompress(payload)), separator=" ")
    print(f"rows {frame.height}, chrom values {frame['CHROM'].unique().to_list()}")
    print(f"fraction of tar fetched for this chromosome: {size / total:.3f}")


if __name__ == "__main__":
    main()
