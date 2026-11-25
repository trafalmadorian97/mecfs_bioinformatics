"""
Adapter Task: convert an LDSC-SEG style directory (a .ldcts index plus chromosome-split LD-score
files) into a prepared S-LDSC annotation set (see prepared_sldsc_reference).

Each .ldcts row is "label<TAB>target_prefix,control_prefix". Every unique prefix is read once,
so a control shared by every row is stored once and referenced by every model.

A prefix may hold several LD columns (the Roadmap control of Corces_ATAC holds 12). As in ldsc's
cell_type_specific, a model's design is the target prefix's columns followed by the control
prefix's columns, and the reported coefficient is the first column's; so the target is the first
column of the target prefix and every other column is a control.
"""

from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_cell_analysis_by_sldsc import (
    LDCTSFileEntry,
    ldcts_raw_to_processed,
)
from mecfs_bio.build_system.task.sldsc.ld_score_files import read_ld_score_sources
from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    LdScoreTable,
    PreparedAnnotationSet,
    SLDSCAnnotationModel,
    write_prepared_annotation_set,
)
from mecfs_bio.build_system.wf.base_wf import WF

PREPARED_SLDSC_SUB_FOLDER = PurePath("prepared_sldsc")


@frozen(slots=True)
class PrepareLdctsAnnotationSetTask(Task):
    """
    source_directory_task: directory containing ldcts_filename and the LD-score files its
        rows reference (paths in the index are relative to this directory).
    ldcts_filename: name of the .ldcts index inside that directory.
    """

    meta: Meta
    source_directory_task: Task
    ldcts_filename: str

    @property
    def deps(self) -> list[Task]:
        return [self.source_directory_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.source_directory_task.asset_id)
        assert isinstance(source, DirectoryAsset)
        entries = _read_ldcts_entries(source.path / self.ldcts_filename)
        prefixes = _unique_prefixes(entries)
        ld = read_ld_score_sources(
            root=source.path, prefixes=prefixes, drop_zero_variance=False
        )
        columns_of = _columns_by_prefix(ld=ld, prefixes=prefixes)
        annotation_set = PreparedAnnotationSet(
            ld=ld, models=tuple(_model_from_entry(e, columns_of) for e in entries)
        )
        out_dir = scratch_dir / "prepared_annotation_set"
        write_prepared_annotation_set(out_dir, annotation_set)
        return DirectoryAsset(out_dir)

    @classmethod
    def create(
        cls, asset_id: str, source_directory_task: Task, ldcts_filename: str
    ) -> "PrepareLdctsAnnotationSetTask":
        source_meta = source_directory_task.meta
        assert isinstance(source_meta, ReferenceDataDirectoryMeta)
        return cls(
            meta=ReferenceDataDirectoryMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PREPARED_SLDSC_SUB_FOLDER,
                id=AssetId(asset_id),
            ),
            source_directory_task=source_directory_task,
            ldcts_filename=ldcts_filename,
        )


def _read_ldcts_entries(path: Path) -> list[LDCTSFileEntry]:
    assert path.is_file(), f"missing .ldcts index {path}"
    entries = ldcts_raw_to_processed(path.read_text().strip())
    assert entries, f"empty .ldcts index {path}"
    return entries


def _unique_prefixes(entries: list[LDCTSFileEntry]) -> list[str]:
    """Every target and control prefix in the index, each listed once, in order of first appearance."""
    prefixes: list[str] = []
    for entry in entries:
        for prefix in (
            str(entry.target_gene_data_path),
            str(entry.control_gene_data_path),
        ):
            if prefix not in prefixes:
                prefixes.append(prefix)
    return prefixes


def _columns_by_prefix(ld: LdScoreTable, prefixes: list[str]) -> dict[str, list[str]]:
    """Column ids of each prefix, in file column order (sources are "prefix:column")."""
    columns_of: dict[str, list[str]] = {prefix: [] for prefix in prefixes}
    for column_id, source in zip(ld.column_ids, ld.sources):
        columns_of[source.rpartition(":")[0]].append(column_id)
    assert all(columns_of.values()), "every cts prefix must hold at least one LD column"
    return columns_of


def _model_from_entry(
    entry: LDCTSFileEntry, columns_of: dict[str, list[str]]
) -> SLDSCAnnotationModel:
    target, *extra_targets = columns_of[str(entry.target_gene_data_path)]
    return SLDSCAnnotationModel(
        label=entry.label,
        target_column=target,
        control_columns=(
            *extra_targets,
            *columns_of[str(entry.control_gene_data_path)],
        ),
    )
