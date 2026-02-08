from pathlib import Path

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
from mecfs_bio.build_system.meta.meta import Meta
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
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.rebuilder.metadata_to_path.base_meta_to_path import (
    MetaToPath,
)


@frozen
class SimpleMetaToPath(MetaToPath):
    root: Path

    def __call__(self, m: Meta) -> Path:
        if isinstance(m, SimpleFileMeta):
            return self.root / "other_files" / m.id
        if isinstance(m, SimpleDirectoryMeta):
            return self.root / "other_files" / m.asset_id
        if isinstance(m, GWASSummaryDataFileMeta):
            pth = self.root / "gwas" / m.trait / m.project / m.sub_dir
            if m.project_path is not None:
                pth = pth / m.project_path
            else:
                f_name = str(m.id)
                if m.extension is not None:
                    f_name += m.extension
                pth = pth / f_name
            return pth
        if isinstance(m, FilteredGWASDataMeta):
            pth = (
                self.root
                / "gwas"
                / m.trait
                / m.project
                / m.sub_dir
                / str(m.id + m.extension)
            )
            return pth
        if isinstance(m, ProcessedGwasDataDirectoryMeta):
            pth = self.root / "gwas" / m.trait / m.project / m.sub_dir / str(m.id)
            return pth
        if isinstance(m, GWASLabSumStatsMeta):
            pth = (
                self.root
                / "gwas"
                / m.trait
                / m.project
                / m.sub_dir
                / (m.asset_id + ".pickle")
            )
            return pth
        if isinstance(m, GWASLabLeadVariantsMeta):
            pth = (
                self.root
                / "gwas"
                / m.trait
                / m.project
                / m.sub_dir
                / str(m.asset_id + ".csv")
            )
            return pth
        if isinstance(m, GWASLabRegionPlotsMeta):
            pth = self.root / "gwas" / m.trait / m.project / m.sub_dir / m.asset_id
            return pth
        if isinstance(m, GWASLabManhattanQQPlotMeta):
            pth = (
                self.root
                / "gwas"
                / m.trait
                / m.project
                / m.sub_dir
                / str(m.asset_id + ".png")
            )
            return pth
        if isinstance(m, ReferenceFileMeta):
            pth = self.root / "reference_data" / m.group / m.sub_group / m.sub_folder
            if m.filename is not None:
                pth = pth / (m.filename + m.extension)
            else:
                pth = pth / str(m.id + m.extension)
            return pth

        if isinstance(m, ReferenceDataDirectoryMeta):
            dirname = m.dirname if m.dirname is not None else m.id
            pth = (
                self.root
                / "reference_data"
                / m.group
                / m.sub_group
                / m.sub_folder
                / dirname
            )
            return pth

        if isinstance(m, ExecutableMeta):
            pth = self.root / "executable" / m.group / m.sub_folder
            if m.filename is not None:
                fname = m.filename
                if m.extension is not None:
                    fname += m.extension
                pth = pth / fname
            else:
                pth = pth / m.id
            return pth

        if isinstance(m, (GWASPlotDirectoryMeta)):
            pth = self.root / "gwas" / m.trait / m.project / m.sub_dir / m.asset_id
            return pth

        if isinstance(m, ( GWASPlotFileMeta)):
            pth = self.root / "gwas" / m.trait / m.project / m.sub_dir / (m.asset_id+m.extension)
            return pth
        if isinstance(m, ResultTableMeta):
            pth = (
                self.root
                / "gwas"
                / m.trait
                / m.project
                / m.sub_dir
                / (m.id + m.extension)
            )
            return pth
        if isinstance(m, ResultDirectoryMeta):
            pth = self.root / "gwas" / m.trait / m.project / m.sub_dir / m.id
            return pth
        raise ValueError(f"Unknown meta {m} of type {type(m)}.")
