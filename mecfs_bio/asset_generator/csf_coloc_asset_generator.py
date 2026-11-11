"""
Asset generator for colocalizing one trait GWAS with the Western et al. 2024 CSF aptamers: the
CSF source (HapMap3 CSF hit table, remote-tabix window fetcher over the GWAS Catalog harmonised
files, Entrez-id gene coordinates from the GRCh38 MAGMA gene table) plugged into the shared
protein coloc generator.

The pair budget is larger than PPP's (CSF_SELECTION: 200 pairs, 20 per locus). PPP's budget
exists because one PPP pair costs a whole Synapse tar download; one CSF pair costs a tabix
index download plus a range read of a few seconds, so 200 pairs take roughly 15-30 minutes.
The per-locus cap is raised because several aptamers often target one gene.
"""

from pathlib import Path

from mecfs_bio.asset_generator.csf_slim_aptamer_asset_generator import (
    CSF_APTAMER_MANIFEST_PATH,
    load_csf_aptamer_files,
)
from mecfs_bio.asset_generator.protein_coloc_asset_generator import (
    ColocTaskGroup,
    ProteomicsColocSource,
    generate_protein_coloc_tasks,
)
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_hit_table import (
    HAPMAP_3_CSF_HIT_TABLE_5E8,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_38_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import TraitModel
from mecfs_bio.build_system.task.protein_coloc.csf_gene_coords_reader import (
    CsfGeneCoordsReader,
)
from mecfs_bio.build_system.task.protein_coloc.csf_window_fetcher import (
    CsfWindowFetcher,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)

CSF_SELECTION = SelectionConfig(max_pairs=200, max_per_locus=20)
CSF_COLOC_CONFIG = ColocRunConfig(selection=CSF_SELECTION)


def generate_csf_coloc_tasks(
    base_name: str,
    trait_task: Task,
    trait_model: TraitModel,
    config: ColocRunConfig = CSF_COLOC_CONFIG,
    h4_threshold: float = 0.8,
    manifest_path: Path = CSF_APTAMER_MANIFEST_PATH,
) -> ColocTaskGroup:
    """
    base_name: prefix of the asset ids ({base_name}_csf_coloc, ...).
    trait_task: the trait GWAS, genome-reference harmonized to GRCh38.
    trait_model: quantitative or case-control, as coloc.abf needs.
    """
    aptamers = load_csf_aptamer_files(manifest_path)
    source = ProteomicsColocSource(
        asset_infix="csf",
        hit_table_task=HAPMAP_3_CSF_HIT_TABLE_5E8,
        gene_coords_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
        gene_coords_reader=CsfGeneCoordsReader(aptamers=aptamers),
        fetcher=CsfWindowFetcher(aptamers=aptamers),
    )
    return generate_protein_coloc_tasks(
        base_name=base_name,
        source=source,
        trait_task=trait_task,
        trait_model=trait_model,
        config=config,
        h4_threshold=h4_threshold,
    )
