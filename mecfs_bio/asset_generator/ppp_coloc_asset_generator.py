"""
Asset generator for colocalizing one trait GWAS with the UKB-PPP proteins.

Builds, for one GRCh38-harmonized trait: the TraitProteinColocTask (raw results and aligned
windows), the docs-facing display table of the pairs with PP.H4 above a threshold, and the
stacked locus plots of those pairs. The shared inputs (protein hit table, PPP manifest, ST3
gene coordinates, GRCh38 genome, allele-frequency panel and gene track) are defaulted here so
each trait's module is a single call.
"""

from pathlib import Path

from attrs import frozen

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
    load_ppp_protein_files,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_hit_table import (
    HAPMAP_3_PPP_HIT_TABLE_5E8,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_38_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
)
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_st3_extracted import (
    SUN_ET_AL_2023_ST3_EXTRACTED,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import TraitModel
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ColocDisplayTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_locus_plot_task import (
    ColocLocusPlotTask,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import (
    PppWindowFetcher,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    PppProteinHitTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
    TraitProteinColocTask,
)


@frozen(slots=True)
class PppColocTaskGroup:
    coloc_task: TraitProteinColocTask
    display_task: ColocDisplayTableTask
    plot_task: ColocLocusPlotTask

    def terminal_tasks(self) -> list[Task]:
        return [self.display_task, self.plot_task]


def generate_ppp_coloc_tasks(
    base_name: str,
    trait_task: Task,
    trait_model: TraitModel,
    config: ColocRunConfig = ColocRunConfig(),
    h4_threshold: float = 0.8,
    hit_table_task: PppProteinHitTableTask = HAPMAP_3_PPP_HIT_TABLE_5E8,
    manifest_path: Path = EUR_DISCOVERY_PPP_MANIFEST_PATH,
) -> PppColocTaskGroup:
    """
    base_name: prefix of the asset ids.
    trait_task: the trait GWAS, genome-reference harmonized to GRCh38.
    trait_model: quantitative or case-control, as coloc.abf needs.
    manifest_path: the PPP cohort whose dense summary statistics are fetched; must match the
        cohort of hit_table_task.
    """
    coloc_task = TraitProteinColocTask.create(
        asset_id=f"{base_name}_ppp_coloc",
        trait_task=trait_task,
        trait_model=trait_model,
        hit_table_task=hit_table_task,
        gene_coords_task=SUN_ET_AL_2023_ST3_EXTRACTED,
        fasta_task=UCSC_HG38_INDEXED_FASTA,
        panel_task=THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
        panel_ancestry="eur",
        fetcher=PppWindowFetcher(proteins=load_ppp_protein_files(manifest_path)),
        config=config,
    )
    display_task = ColocDisplayTableTask.create(
        asset_id=f"{base_name}_ppp_coloc_display",
        coloc_task=coloc_task,
        h4_threshold=h4_threshold,
    )
    plot_task = ColocLocusPlotTask.create(
        asset_id=f"{base_name}_ppp_coloc_plots",
        coloc_task=coloc_task,
        display_task=display_task,
        gene_info_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
    )
    return PppColocTaskGroup(
        coloc_task=coloc_task, display_task=display_task, plot_task=plot_task
    )
