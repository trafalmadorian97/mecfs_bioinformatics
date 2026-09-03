"""Build one sorted, allele-bearing annotation parquet from the baseline-LF
per-chromosome .annot.parquet members.

Consumes the directory of baselineLF2.2.UKB.<chr>.annot.parquet members produced
by StreamExtractAnnotationParquetsTask (each: CHR, SNP, BP, A1, A2, + 187
annotation columns; no CM). Casts the annotation columns to float32, collapses
the rare unordered-allele ordering duplicates (same position with A1/A2 swapped
and identical annotations), and streams them, in (CHR, BP) order, into a single
sorted parquet. The result is keyed allele-aware: it is unique on
(CHR, BP, unordered-allele-key). Downstream tasks (ridge weights, explainability)
join to gwas variants on (CHR, BP, unordered-allele-key), so each allele of a
multiallelic site carries its own annotations.
"""

import re
from pathlib import Path

import polars as pl
import pyarrow.parquet as pq
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key
from mecfs_bio.build_system.wf.base_wf import WF

# The annotation source now carries alleles, so the position key includes A1/A2.
ANNOT_KEY_COLUMNS: list[str] = ["CHR", "BP", "SNP", "A1", "A2"]
_CHR_COL = "CHR"
_BP_COL = "BP"
_A1_COL = "A1"
_A2_COL = "A2"
_ALLELE_KEY_COL = "allele_key"
_ANNOT_MEMBER_RE = re.compile(r"baselineLF2\.2\.UKB\.(\d+)\.annot\.parquet$")


@frozen
class BuildBaselineLFAnnotationParquetTask(Task):
    """Build a single sorted allele-bearing annotation parquet from the members."""

    meta: Meta
    annot_members_task: Task

    @property
    def deps(self) -> list["Task"]:
        return [self.annot_members_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        members_asset = fetch(self.annot_members_task.asset_id)
        assert isinstance(members_asset, DirectoryAsset)
        members = _list_annot_members(members_asset.path)

        # Process one chromosome at a time and append it to a single parquet via
        # a streaming pyarrow writer. A global concat/sort/unique over all ~19.5M
        # rows OOMs the 16GB budget, so we never hold more than one deduped
        # ~1.5M-row chromosome in memory. Each member is already BP-sorted, so
        # writing chromosomes in order yields a globally (CHR, BP)-sorted result.
        out_path = scratch_dir / str(self.meta.asset_id)
        writer: pq.ParquetWriter | None = None
        try:
            for chrom in sorted(members):
                table = _dedup_one_chromosome(members[chrom]).to_arrow()
                if writer is None:
                    writer = pq.ParquetWriter(out_path, table.schema)
                writer.write_table(table)
        finally:
            if writer is not None:
                writer.close()
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, annot_members_task: Task
    ) -> "BuildBaselineLFAnnotationParquetTask":
        # Derive the output file meta from the members directory dependency's
        # meta (reuse group/sub_group/sub_folder); only id/extension/read_spec
        # differ.
        source_meta = annot_members_task.meta
        if not isinstance(source_meta, ReferenceDataDirectoryMeta):
            raise ValueError(f"Unknown meta for annotation members task: {source_meta}")
        meta = ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta=meta, annot_members_task=annot_members_task)


def _list_annot_members(members_dir: Path) -> dict[int, Path]:
    """Map chromosome -> member path for every baselineLF *.annot.parquet file."""
    members: dict[int, Path] = {}
    for path in members_dir.iterdir():
        match = _ANNOT_MEMBER_RE.search(path.name)
        if match:
            members[int(match.group(1))] = path
    assert members, f"no *.annot.parquet members found in {members_dir}"
    return members


def _dedup_one_chromosome(member_path: Path) -> pl.DataFrame:
    """Scan one chromosome's member, cast annotations to float32, and collapse its
    unordered-allele ordering duplicates (identical annotations), returning an
    eager frame. Done per chromosome (bounded memory) since such duplicates share
    a (CHR, BP); genuine multiallelic sites have distinct allele keys and survive.

    Collapsing on the position key alone relies on all rows sharing a
    (CHR, BP, unordered-allele-key) carrying identical annotations. We enforce
    that cheaply instead of assuming it: dedup on the key AND every annotation
    column, so ordering duplicates that agree collapse, and any key that still
    appears more than once must disagree on some annotation -- which we reject.
    """
    lazy = pl.scan_parquet(member_path)
    schema = lazy.collect_schema()
    annot_cols = [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]
    key_cols = [_CHR_COL, _BP_COL, _ALLELE_KEY_COL]
    deduped = (
        lazy.with_columns([pl.col(c).cast(pl.Float32) for c in annot_cols])
        .with_columns(unordered_allele_key(_A1_COL, _A2_COL).alias(_ALLELE_KEY_COL))
        .unique(subset=[*key_cols, *annot_cols], keep="first")
        .collect()
    )
    conflicting_keys = deduped.height - deduped.n_unique(subset=key_cols)
    assert conflicting_keys == 0, (
        f"{member_path.name}: {conflicting_keys} (CHR, BP, unordered-allele-key) "
        "group(s) carry differing annotations, violating the dedup assumption"
    )
    return (
        deduped.drop(_ALLELE_KEY_COL)
        # unique() may reorder; restore per-chromosome BP order so the streamed
        # concatenation of chromosomes is globally (CHR, BP)-sorted.
        .sort(_BP_COL)
    )
