"""
Asset generator for colocalizing one trait GWAS with the deCODE SomaScan assays (Eldjarn et al.
2023): the deCODE source of one group (HapMap3 hit table, window reads of the remote files
through each protein's stored seek index, Entrez-id gene coordinates from the GRCh38 MAGMA gene
table) plugged into the shared protein coloc generator.

One deCODE pair costs a HEAD request and a few range reads of the protein's remote file, so the
pair budget is a choice of the run, not a cost limit.
"""

from mecfs_bio.asset_generator.decode_somascan_protein_asset_generator import (
    DecodeSomascanTaskCollection,
)
from mecfs_bio.asset_generator.protein_coloc_asset_generator import (
    ColocTaskGroup,
    ProteomicsColocSource,
    ProteomicsSourceInfix,
    generate_protein_coloc_tasks,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_database_protein_files import (
    HAPMAP_3_DECODE_RAW_DATABASE,
    HAPMAP_3_DECODE_SMP_DATABASE,
)
from mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_hit_table import (
    HAPMAP_3_DECODE_RAW_HIT_TABLE_5E8,
    HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8,
)
from mecfs_bio.assets.gwas.decode_somascan.raw.decode_variant_info import (
    DECODE_VARIANT_INFO,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_38_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import TraitModel
from mecfs_bio.build_system.task.protein_coloc.decode_gene_coords_reader import (
    DecodeGeneCoordsReader,
)
from mecfs_bio.build_system.task.protein_coloc.decode_window_reader import (
    DECODE_MAX_PARALLEL_FETCHES,
    DecodeProteinEntry,
    DecodeWindowReader,
)
from mecfs_bio.build_system.task.protein_coloc.dense_window_fetch import (
    ParallelWindowFetcher,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    ProteinHitTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)
from mecfs_bio.constants.decode_somascan_constants import DecodeSomascanGroup

_DATABASES: dict[DecodeSomascanGroup, DecodeSomascanTaskCollection] = {
    "smp": HAPMAP_3_DECODE_SMP_DATABASE,
    "raw": HAPMAP_3_DECODE_RAW_DATABASE,
}
_INFIXES: dict[DecodeSomascanGroup, ProteomicsSourceInfix] = {
    "smp": "decode_somascan_smp",
    "raw": "decode_somascan_raw",
}
_HIT_TABLES: dict[DecodeSomascanGroup, ProteinHitTableTask] = {
    "smp": HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8,
    "raw": HAPMAP_3_DECODE_RAW_HIT_TABLE_5E8,
}


def generate_decode_somascan_coloc_tasks(
    base_name: str,
    trait_task: Task,
    trait_model: TraitModel,
    config: ColocRunConfig,
    group: DecodeSomascanGroup = "smp",
    h4_threshold: float = 0.8,
) -> ColocTaskGroup:
    """
    base_name: prefix of the asset ids ({base_name}_decode_somascan_{group}_coloc, ...).
    trait_task: the trait GWAS, genome-reference harmonized to GRCh38.
    trait_model: quantitative or case-control, as coloc.abf needs.
    group: SMP-normalized ("smp", the default) or non-normalized ("raw") assays.
    """
    protein_tasks = _DATABASES[group].protein_tasks
    source = ProteomicsColocSource(
        asset_infix=_INFIXES[group],
        hit_table_task=_HIT_TABLES[group],
        gene_coords_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
        gene_coords_reader=DecodeGeneCoordsReader(
            proteins=tuple(task.protein for task in protein_tasks)
        ),
        fetcher=ParallelWindowFetcher(
            reader=DecodeWindowReader(
                proteins=tuple(DecodeProteinEntry(task) for task in protein_tasks),
                variant_info_task=DECODE_VARIANT_INFO,
            ),
            max_parallel=DECODE_MAX_PARALLEL_FETCHES,
        ),
    )
    return generate_protein_coloc_tasks(
        base_name=base_name,
        source=source,
        trait_task=trait_task,
        trait_model=trait_model,
        config=config,
        h4_threshold=h4_threshold,
    )
