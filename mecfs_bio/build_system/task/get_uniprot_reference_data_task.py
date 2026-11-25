"""
Task to download reference data about proteins from UniProt
"""

import functools
from pathlib import Path, PurePath
from typing import Self

import requests
from attrs import frozen
from UniProtMapper import ProtKB
from UniProtMapper.uniprotkb_fields import organism_name, reviewed

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.retry.call_with_retries import call_with_retries

protkb = ProtKB()

DEFAULT_FIELDS = [
    "accession",
    "id",
    "gene_names",
    "protein_name",
    "organism_name",
    "organism_id",
    "go_id",
    "go_p",
    "go_c",
    "go_f",
    "cc_subcellular_location",
    "cc_function",
    "xref_ensembl",
    "xref_reactome",
]

# UniProtMapper's urllib3 Retry covers connection failures and error statuses,
# but not a response body truncated mid-stream, which surfaces as
# ChunkedEncodingError while paging. A retry restarts the full paged query.
_TRANSIENT_UNIPROT_ERRORS: tuple[type[Exception], ...] = (
    requests.exceptions.ChunkedEncodingError,
    requests.exceptions.ConnectionError,
)


@frozen(slots=True)
class GetUniProtReferenceDataTask(Task):
    """
    Task to download reference data about proteins from UniProt
    """

    meta: Meta
    field_list: list[str] = (
        DEFAULT_FIELDS  # see: https://david-araripe.github.io/UniProtMapper/stable/field_reference.html
    )

    @property
    def deps(self) -> list["Task"]:
        return []

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        query = reviewed(True) & organism_name("human")  # type: ignore[operator] # ty: ignore[unsupported-operator]
        result = call_with_retries(
            functools.partial(protkb.get, query, fields=self.field_list),
            retry_on=_TRANSIENT_UNIPROT_ERRORS,
        )
        out_path = scratch_dir / "uniprot.parquet"
        result.to_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls, asset_id: str) -> Self:
        meta = ReferenceFileMeta(
            group="protein_lookup",
            sub_group="uniprot",
            sub_folder=PurePath("raw"),
            id=AssetId(asset_id),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta)
