from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.executable_meta import ExecutableMeta
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_lead_variants_meta import (
    GWASLabLeadVariantsMeta,
)
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_manhattan_plot_meta import (
    GWASLabManhattanQQPlotMeta,
)
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_region_plots_meta import (
    GWASLabRegionPlotsMeta,
)
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.markdown_file_meta import MarkdownFileMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.plot_meta import GWASPlotDirectoryMeta
from mecfs_bio.build_system.meta.processed_gwas_data_directory_meta import (
    ProcessedGwasDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.result_archive_meta import ResultArchiveMeta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)

_GWAS = PurePath("gwas")
_REFERENCE_DATA = PurePath("reference_data")
_OTHER_FILES = PurePath("other_files")
_EXECUTABLE = PurePath("executable")


def simple_meta_to_relative_path(m: Meta) -> PurePath:
    """
    Compute where an asset lives inside the asset store, as a path relative to the store
    root.

    The layout is deliberately root-free so that one logical asset store can be spread
    across several filesystems: the same relative path is simply appended to whichever
    root the caller selects.  See the remapping_meta_to_path module for the machinery
    that exploits this.
    """
    if isinstance(m, SimpleFileMeta):
        return _OTHER_FILES / m.id
    if isinstance(m, SimpleDirectoryMeta):
        return _OTHER_FILES / m.asset_id
    if isinstance(m, GWASSummaryDataFileMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir
        if m.project_path is not None:
            pth = pth / m.project_path
        else:
            f_name = str(m.id)
            if m.extension is not None:
                f_name += m.extension
            pth = pth / f_name
        return pth
    if isinstance(m, FilteredGWASDataMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / str(m.id + m.extension)
        return pth
    if isinstance(m, MarkdownFileMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / str(m.id + ".md")
        return pth
    if isinstance(m, ProcessedGwasDataDirectoryMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / str(m.id)
        return pth
    if isinstance(m, GWASLabSumStatsMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / (m.asset_id + ".pickle")
        return pth
    if isinstance(m, GWASLabLeadVariantsMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / str(m.asset_id + ".csv")
        return pth
    if isinstance(m, GWASLabRegionPlotsMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / m.asset_id
        return pth
    if isinstance(m, GWASLabManhattanQQPlotMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / str(m.asset_id + ".png")
        return pth
    if isinstance(m, ReferenceFileMeta):
        pth = _REFERENCE_DATA / m.group / m.sub_group / m.sub_folder
        if m.filename is not None:
            pth = pth / (m.filename + m.extension)
        else:
            pth = pth / str(m.id + m.extension)
        return pth

    if isinstance(m, ReferenceDataDirectoryMeta):
        dirname = m.dirname if m.dirname is not None else m.id
        pth = _REFERENCE_DATA / m.group / m.sub_group / m.sub_folder / dirname
        return pth

    if isinstance(m, ExecutableMeta):
        pth = _EXECUTABLE / m.group / m.sub_folder
        if m.filename is not None:
            fname = m.filename
            if m.extension is not None:
                fname += m.extension
            pth = pth / fname
        else:
            pth = pth / m.id
        return pth

    if isinstance(m, (GWASPlotDirectoryMeta)):
        pth = _GWAS / m.trait / m.project / m.sub_dir / m.asset_id
        return pth

    if isinstance(m, (GWASPlotFileMeta)):
        pth = _GWAS / m.trait / m.project / m.sub_dir / (m.asset_id + m.extension)
        return pth
    if isinstance(m, ResultTableMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / (m.id + m.extension)
        return pth
    if isinstance(m, ResultDirectoryMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / m.id
        return pth
    if isinstance(m, ResultArchiveMeta):
        pth = _GWAS / m.trait / m.project / m.sub_dir / (m.id + m.extension)
        return pth
    raise ValueError(f"Unknown meta {m} of type {type(m)}.")


@frozen
class SimpleMetaToPath(MetaToPath):
    """
    Place every asset under a single store root, using the standard relative layout.
    """

    root: Path

    def __call__(self, m: Meta) -> Path:
        return self.root / simple_meta_to_relative_path(m)
