"""
A gnomAD release's allele-frequency panel: the per-chromosome tables concatenated in gwaslab
chromosome order into one parquet sorted by (CHR, POS), so the harmonizer's per-chromosome
filter prunes row groups.

The per-chromosome Tasks are injected by the caller (normally
generate_gnomad_allele_frequency_panel_tasks) and are ordinary dependencies, kept in the asset
store under their own sub_folder. They cannot be deleted to save space, because the build system materializes
every transitive dependency of a target; a path_remap rule can move them to another disk.
Do not wrap this Task in DiscardDepsWrapper: the multi-hour build would become
all-or-nothing again and the FASTA would be rebuilt in a temporary store.
"""

from pathlib import Path, PurePath
from typing import Sequence

import pyarrow
import pyarrow.parquet
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    open_parquet_writer,
    parquet_encoding,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GNOMAD_GROUP,
    GnomadRelease,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    fetch_file_path,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import panel_an_col
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL

GNOMAD_PANEL_FILENAME = "gnomad_allele_frequencies"


@frozen(slots=True)
class GnomadAlleleFrequencyPanelTask(Task):
    meta: HarmonizableReferenceTableMeta
    release: GnomadRelease
    part_tasks: tuple[GnomadChromosomeAlleleFrequencyTask, ...]

    def __attrs_post_init__(self) -> None:
        # Exactly the release's chromosomes, in release order: concatenation relies on it
        # for a complete panel sorted by (CHR, POS).
        part_chromosomes = tuple(part.chrom for part in self.part_tasks)
        assert part_chromosomes == self.release.chromosomes, (
            f"{self.release.name}: parts cover chromosomes {part_chromosomes}, "
            f"expected {self.release.chromosomes} in that order"
        )
        foreign = [
            part.asset_id for part in self.part_tasks if part.release != self.release
        ]
        assert not foreign, (
            f"parts from a release other than {self.release.name}: {foreign}"
        )

    @property
    def deps(self) -> list[Task]:
        return list(self.part_tasks)

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        part_paths = [fetch_file_path(fetch, part) for part in self.part_tasks]
        schema = identical_part_schema(part_paths)
        out_path = scratch_dir / (GNOMAD_PANEL_FILENAME + ".parquet")
        concatenate_parts(
            part_paths,
            out_path,
            schema,
            byte_stream_split_columns=[
                panel_an_col(group) for group in self.release.groups
            ],
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        release: GnomadRelease,
        part_tasks: Sequence[GnomadChromosomeAlleleFrequencyTask],
    ) -> "GnomadAlleleFrequencyPanelTask":
        return cls(
            meta=HarmonizableReferenceTableMeta(
                group=GNOMAD_GROUP,
                sub_group=release.name,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename=GNOMAD_PANEL_FILENAME,
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=release.build,
                    ref_allele_col=PANEL_REF_COL,
                    pos_col=GWASLAB_POS_COL,
                ),
                allele_frequency_columns=PanelAlleleFrequencyColumns.prefixed(
                    release.groups
                ),
            ),
            release=release,
            part_tasks=tuple(part_tasks),
        )


def identical_part_schema(part_paths: list[Path]) -> pyarrow.Schema:
    """The parts' shared arrow schema (metadata ignored); fails if any part differs."""
    schemas = [
        pyarrow.parquet.read_schema(path).remove_metadata() for path in part_paths
    ]
    assert schemas, "no per-chromosome parts"
    differing = [
        str(path)
        for path, schema in zip(part_paths, schemas)
        if not schema.equals(schemas[0])
    ]
    assert not differing, (
        f"parts whose schema differs from {part_paths[0]}: {differing}"
    )
    return schemas[0]


def concatenate_parts(
    part_paths: list[Path],
    out_path: Path,
    schema: pyarrow.Schema,
    byte_stream_split_columns: list[str],
) -> None:
    """Copy every row group of every part, in order, into one parquet."""
    encoding = parquet_encoding(
        schema.names,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=byte_stream_split_columns,
    )
    with open_parquet_writer(out_path, schema, encoding) as writer:
        for path in part_paths:
            part = pyarrow.parquet.ParquetFile(path)
            for index in range(part.num_row_groups):
                writer.write_table(
                    part.read_row_group(index).replace_schema_metadata(None)
                )
