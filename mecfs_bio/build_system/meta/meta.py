"""
Meta objects are uniquely identifying metadata describing either an asset that currently exists, or an asset that can be created by a build system.
"""

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

Meta = (
    SimpleFileMeta
    | SimpleDirectoryMeta
    | GWASSummaryDataFileMeta
    | FilteredGWASDataMeta
    | GWASLabSumStatsMeta
    | GWASLabLeadVariantsMeta
    | GWASLabRegionPlotsMeta
    | ReferenceFileMeta
    | GWASLabManhattanQQPlotMeta
    | ReferenceDataDirectoryMeta
    | ExecutableMeta
    | ProcessedGwasDataDirectoryMeta
    | GWASPlotDirectoryMeta
    | ResultTableMeta
    | ResultDirectoryMeta
)
