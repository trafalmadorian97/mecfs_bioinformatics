"""Build one sorted annotation parquet from the baseline-LF 2.2.UKB tarball.

Extracts the per-chromosome baselineLF2.2.UKB.<chr>.annot.gz members (LDSC text:
CHR, BP, SNP, CM + annotation columns), casts the annotation columns to float32,
dedups multiallelic rsid duplicates, and streams them, in chromosome order, into
a single (CHR, BP)-sorted parquet. The single sorted file gives the ridge weights
task a cheap scan and the downstream explainability tasks predicate-pushdown
per-locus lookups.
"""

import re
from pathlib import Path

import narwhals
import polars as pl
import structlog
from attrs import frozen

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
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    write_df_according_to_format,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

ANNOT_KEY_COLUMNS: list[str] = ["CHR", "BP", "SNP", "CM"]
_SNP_COL = "SNP"
_CHR_COL = "CHR"
_BP_COL = "BP"
_ANNOT_MEMBER_RE = re.compile(r"baselineLF2\.2\.UKB\.(\d+)\.annot\.gz$")


@frozen
class BuildBaselineLFAnnotationParquetTask(Task):
    """Build a single sorted, deduped annotation parquet from the baseline-LF tarball."""

    meta: Meta
    tarball_task: Task

    @property
    def deps(self) -> list["Task"]:
        return [self.tarball_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        tarball_asset = fetch(self.tarball_task.asset_id)
        assert isinstance(tarball_asset, FileAsset)
        members = _list_annot_members(tarball_asset.path)
        extract_dir = scratch_dir / "annot_gz"
        extract_dir.mkdir(parents=True, exist_ok=True)
        _extract_members(tarball_asset.path, members, extract_dir)

        lazy_frames = []
        for chrom in sorted(members):
            member_path = extract_dir / members[chrom]
            lazy_frames.append(_scan_one_chromosome(member_path))
        combined = pl.concat(lazy_frames, how="vertical").sort([_CHR_COL, _BP_COL])

        out_path = scratch_dir / str(self.meta.asset_id)
        write_df_according_to_format(
            df=narwhals.from_native(combined).lazy(),
            out_path=out_path,
            out_format=ParquetOutFormat(),
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, tarball_task: Task
    ) -> "BuildBaselineLFAnnotationParquetTask":
        # Derive the output meta from the tarball dependency's meta (reuse
        # group/sub_group/sub_folder), the CompressedCSVToParquetTask.create
        # pattern. Only the id, extension and read_spec differ.
        source_meta = tarball_task.meta
        if not isinstance(source_meta, ReferenceFileMeta):
            raise ValueError(f"Unknown meta for tarball task: {source_meta}")
        meta = ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta=meta, tarball_task=tarball_task)


def _list_annot_members(tarball: Path) -> dict[int, str]:
    """Map chromosome -> archive member path for every *.annot.gz member."""
    listing = execute_command(["tar", "-tzf", str(tarball)])
    members: dict[int, str] = {}
    for line in listing.splitlines():
        name = line.strip()
        match = _ANNOT_MEMBER_RE.search(name)
        if match:
            members[int(match.group(1))] = name
    assert members, f"no *.annot.gz members found in {tarball}"
    return members


def _extract_members(tarball: Path, members: dict[int, str], dest: Path) -> None:
    execute_command(["tar", "-xzf", str(tarball), "-C", str(dest), *members.values()])


def _scan_one_chromosome(member_path: Path) -> pl.LazyFrame:
    lazy = pl.scan_csv(member_path, separator="\t", infer_schema_length=None)
    schema = lazy.collect_schema()
    annot_cols = [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]
    return lazy.with_columns([pl.col(c).cast(pl.Float32) for c in annot_cols]).unique(
        subset=_SNP_COL, keep="first", maintain_order=True
    )
