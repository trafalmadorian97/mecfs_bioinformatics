"""
Asset generator for colocalizing one trait GWAS with the UKB-PPP proteins: the PPP source
(protein hit table, PPP manifest fetcher, Sun et al. 2023 ST3 gene coordinates) plugged into
the shared protein coloc generator, so each trait's module is a single call.
"""

from pathlib import Path

from mecfs_bio.asset_generator.protein_coloc_asset_generator import (
    ColocTaskGroup,
    ProteomicsColocSource,
    generate_protein_coloc_tasks,
)
from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_hit_table import (
    HAPMAP_3_PPP_HIT_TABLE_5E8,
)
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_st3_extracted import (
    SUN_ET_AL_2023_ST3_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import TraitModel
from mecfs_bio.build_system.task.protein_coloc.gene_coords_reader import (
    St3GeneCoordsReader,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import (
    PppWindowFetcher,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    ProteinHitTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
)


def generate_ppp_coloc_tasks(
    base_name: str,
    trait_task: Task,
    trait_model: TraitModel,
    config: ColocRunConfig = ColocRunConfig(),
    h4_threshold: float = 0.8,
    hit_table_task: ProteinHitTableTask = HAPMAP_3_PPP_HIT_TABLE_5E8,
    manifest_path: Path = EUR_DISCOVERY_PPP_MANIFEST_PATH,
) -> ColocTaskGroup:
    """
    base_name: prefix of the asset ids.
    trait_task: the trait GWAS, genome-reference harmonized to GRCh38.
    trait_model: quantitative or case-control, as coloc.abf needs.
    manifest_path: the PPP cohort whose dense summary statistics are fetched; must match the
        cohort of hit_table_task.
    """
    source = ProteomicsColocSource(
        asset_infix="ppp",
        hit_table_task=hit_table_task,
        gene_coords_task=SUN_ET_AL_2023_ST3_EXTRACTED,
        gene_coords_reader=St3GeneCoordsReader(),
        fetcher=PppWindowFetcher(proteins=load_ppp_protein_files(manifest_path)),
    )
    return generate_protein_coloc_tasks(
        base_name=base_name,
        source=source,
        trait_task=trait_task,
        trait_model=trait_model,
        config=config,
        h4_threshold=h4_threshold,
    )
