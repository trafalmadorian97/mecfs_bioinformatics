"""
Recover the per-protein GWAS sample size (N) for every UKB-PPP protein cheaply.

The slim per-protein parquet files store only beta/se, having discarded N; and no
per-protein N is published in the Sun et al. 2023 supplementary tables. But N is a single
constant across all variants of a protein's regenie output, and each protein's raw upload
is a tar of per-chromosome gzipped regenie files (the first file's compressed data starts at
byte 1024). Synapse serves those tars from S3 with HTTP Range support, so we can pull just
the first ~64 KB, parse the tar header, gzip-decompress the prefix, and read N from the
first data row -- ~64 KB per protein instead of a ~550 MB download.

This produces one small table with one row per protein: (oid, gene, synapse_id, n).
"""

from __future__ import annotations

import zlib
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen
from tqdm import tqdm

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import GeneratingTask, Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_N_GENE_COL,
    PPP_N_OID_COL,
    PPP_N_SAMPLE_SIZE_COL,
    PPP_N_SYNID_COL,
)
from mecfs_bio.constants.regenie_constants import REGENIE_N_COL

logger = structlog.get_logger()


@frozen
class PppProteinRef:
    """Identity of one UKB-PPP protein needed to fetch its sample size: the OID primary
    key, gene symbol, and Synapse entity id of the raw tar. The manifest-reading that
    produces these lives in the asset_generator layer, keeping this task manifest-agnostic."""

    oid: str
    gene: str
    synid: str


# 64 KiB is far more than enough: it spans the tar's 512-byte directory stub + the first
# gzip member's header and decompresses to ~1800 regenie rows, while a single Range request
# stays tiny.
DEFAULT_HEAD_BYTES = 64 * 1024

_TAR_BLOCK = 512
_TAR_TYPEFLAG_OFFSET = 156
_TAR_SIZE_OFFSET = 124
_TAR_SIZE_LEN = 12
_TAR_REGULAR_TYPEFLAGS = (b"0", b"\x00")
# gzip streams need wbits = 16 + MAX_WBITS in zlib.
_GZIP_WBITS = 16 + zlib.MAX_WBITS


def _first_regular_member_data_offset(head: bytes) -> int:
    """Offset of the first regular-file member's DATA within a tar buffer (skipping the
    leading directory stub). Tar headers are 512-byte blocks; each member's data follows
    its header, padded up to the next 512-byte boundary."""
    offset = 0
    while offset + _TAR_BLOCK <= len(head):
        header = head[offset : offset + _TAR_BLOCK]
        if header == b"\x00" * _TAR_BLOCK:
            break
        typeflag = header[_TAR_TYPEFLAG_OFFSET : _TAR_TYPEFLAG_OFFSET + 1]
        size_field = header[_TAR_SIZE_OFFSET : _TAR_SIZE_OFFSET + _TAR_SIZE_LEN].rstrip(
            b"\x00 "
        )
        size = int(
            size_field or b"0", 8
        )  # size field is a string representing an octal (base 8) number.  int(...,8) converts this base 8 number to an int
        data_offset = offset + _TAR_BLOCK
        if typeflag in _TAR_REGULAR_TYPEFLAGS and size > 0:
            return data_offset
        offset = data_offset + ((size + _TAR_BLOCK - 1) // _TAR_BLOCK) * _TAR_BLOCK
    raise ValueError("no regular-file member found in tar head buffer")


def extract_regenie_n_from_tar_head(head: bytes) -> int:
    """Read the (constant) regenie sample size N from the first ~KB of a protein tar.

    head: the first bytes of a UKB-PPP protein tar (a directory stub followed by
    per-chromosome gzipped regenie files). Enough bytes must be present to decode the
    tar header plus the first two lines of the first gzip member.
    """
    data_offset = _first_regular_member_data_offset(head)
    decompressed = zlib.decompressobj(_GZIP_WBITS).decompress(head[data_offset:])
    lines = decompressed.decode("utf-8", "replace").splitlines()
    assert len(lines) >= 2, "not enough decompressed regenie rows to read N"
    columns = lines[0].split()
    assert REGENIE_N_COL in columns, (
        f"regenie header missing {REGENIE_N_COL}: {columns}"
    )
    n_index = columns.index(REGENIE_N_COL)
    first_row = lines[1].split()
    return int(float(first_row[n_index]))


@frozen
class PppProteinSampleSizeTask(GeneratingTask):
    """One table of per-protein GWAS sample sizes, recovered by ranged reads from Synapse.

    protein_refs: (oid, gene, synid) for every protein, supplied by the manifest-reading
        generator in the asset_generator layer.
    head_bytes: how many leading bytes of each protein tar to fetch.
    """

    meta: Meta
    protein_refs: tuple[PppProteinRef, ...]
    head_bytes: int

    @property
    def deps(self) -> list[Task]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        records = []
        for ref in tqdm(self.protein_refs):
            head = wf.fetch_synapse_file_head(ref.synid, self.head_bytes)
            n = extract_regenie_n_from_tar_head(head)
            records.append(
                {
                    PPP_N_OID_COL: ref.oid,
                    PPP_N_GENE_COL: ref.gene,
                    PPP_N_SYNID_COL: ref.synid,
                    PPP_N_SAMPLE_SIZE_COL: n,
                }
            )
        table = pl.DataFrame(
            records,
            schema={
                PPP_N_OID_COL: pl.String,
                PPP_N_GENE_COL: pl.String,
                PPP_N_SYNID_COL: pl.String,
                PPP_N_SAMPLE_SIZE_COL: pl.Int64,
            },
        )
        logger.info("recovered ppp per-protein sample sizes", n_proteins=table.height)
        out_path = scratch_dir / f"{self.meta.asset_id}.parquet"
        table.write_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        protein_refs: tuple[PppProteinRef, ...],
        head_bytes: int = DEFAULT_HEAD_BYTES,
    ) -> PppProteinSampleSizeTask:
        meta = ResultTableMeta(
            id=asset_id,
            trait="ukbb_ppp",
            project="ppp_sample_size",
            sub_dir=PurePath("analysis"),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta=meta, protein_refs=protein_refs, head_bytes=head_bytes)
