"""
Time the phases of a ranged PPP window fetch: header walk vs chromosome-file read.

Compares the header walk with a fresh connection per request against a keep-alive
requests.Session, and times reading one large chromosome file over one connection vs several
parallel range requests. Then times the production reader (SharedClientSynapseDownloader,
which keeps one session per thread) walking the headers of three proteins in turn, to show
connections carrying over from one protein to the next.

Run: pixi r python experiments/claude/protein_coloc/time_ppp_ranged_read_phases.py 2>&1 | tee experiments/claude/protein_coloc/time_ppp_ranged_read_phases.log
"""

import time
from concurrent.futures import ThreadPoolExecutor

import requests
import synapseclient
from attrs import frozen

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)
from mecfs_bio.build_system.wf.synapse_downloader import SharedClientSynapseDownloader
from mecfs_bio.util.file_io.ranged_tar import TarMemberExtent, list_tar_members

_N_PARALLEL = (1, 8, 4, 1, 4, 8, 8, 1, 4)


@frozen(slots=True)
class _UrlReader:
    url: str
    session: requests.Session | None

    def read_range(self, start: int, length: int) -> bytes:
        getter = self.session.get if self.session is not None else requests.get
        response = getter(
            self.url,
            headers={"Range": f"bytes={start}-{start + length - 1}"},
            timeout=120,
        )
        assert response.status_code == 206, response.status_code
        return response.content


def _parallel_read(reader: _UrlReader, member: TarMemberExtent, n_parts: int) -> int:
    part = -(-member.size // n_parts)
    starts = range(member.data_offset, member.data_offset + member.size, part)
    end = member.data_offset + member.size
    with ThreadPoolExecutor(max_workers=n_parts) as pool:
        chunks = pool.map(lambda s: reader.read_range(s, min(part, end - s)), starts)
        return sum(len(c) for c in chunks)


def main() -> None:
    protein = load_ppp_protein_files(EUR_DISCOVERY_PPP_MANIFEST_PATH)[0]
    client = synapseclient.login()
    entity = client.get(protein.synid, downloadFile=False)
    url = client._getFileHandleDownload(entity.dataFileHandleId, protein.synid)[
        "preSignedURL"
    ]
    for label, session in (("no session", None), ("keep-alive session", requests.Session())):
        start = time.monotonic()
        members = list_tar_members(_UrlReader(url=url, session=session))
        print(f"header walk, {label}: {time.monotonic() - start:.2f}s "
              f"({len(members)} members)")
    biggest = max(members, key=lambda m: m.size)
    print(f"largest member {biggest.size / 1e6:.1f} MB: {biggest.name}")
    for n_parts in _N_PARALLEL:
        reader = _UrlReader(url=url, session=None)
        start = time.monotonic()
        n_bytes = _parallel_read(reader, biggest, n_parts)
        assert n_bytes == biggest.size
        elapsed = time.monotonic() - start
        print(f"read with {n_parts} connection(s): {elapsed:.2f}s "
              f"({biggest.size / 1e6 / elapsed:.1f} MB/s)")

    _time_production_reader()


def _time_production_reader() -> None:
    downloader = SharedClientSynapseDownloader()
    for protein in load_ppp_protein_files(EUR_DISCOVERY_PPP_MANIFEST_PATH)[1:4]:
        start = time.monotonic()
        reader = downloader.open_range_reader(protein.synid)
        opened = time.monotonic()
        members = list_tar_members(reader)
        print(f"production reader, {protein.oid}: open {opened - start:.2f}s, "
              f"header walk {time.monotonic() - opened:.2f}s ({len(members)} members)")


if __name__ == "__main__":
    main()
