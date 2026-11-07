"""
Asset generator for colocalizing one trait GWAS with the proteins of one proteomics GWAS
(UKB-PPP plasma, Western et al. CSF).

Builds, for one GRCh38-harmonized trait: the TraitProteinColocTask (raw results and aligned
windows) and the docs-facing display table of the pairs with PP.H4 above a threshold. Stacked
locus plots are not built per run: generate_coloc_locus_plot_task makes them for hand-picked
pairs of a run. The proteomics source supplies the hit table, the dense
window fetcher and the gene coordinates; the GRCh38 genome, allele-frequency panel and gene
track are shared by every source.
"""

from typing import Literal

from attrs import frozen

from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_38_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import TraitModel
from mecfs_bio.build_system.task.protein_coloc.coloc_display_table_task import (
    ColocDisplayTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_locus_plot_task import (
    ColocLocusPlotTask,
    ColocPair,
)
from mecfs_bio.build_system.task.protein_coloc.gene_coords_reader import (
    GeneCoordsReader,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import (
    DenseWindowFetcher,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    ProteinHitTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
    TraitProteinColocTask,
)

ProteomicsSourceInfix = Literal["ppp", "csf"]


@frozen(slots=True)
class ProteomicsColocSource:
    """
    asset_infix: names the source in the asset ids ({base}_{infix}_coloc, ...).
    hit_table_task: the source's protein hit-segment table.
    gene_coords_task, gene_coords_reader: the GRCh38 gene-coordinate table and its reader.
    fetcher: dense window source of the proteomics GWAS.
    """

    asset_infix: ProteomicsSourceInfix
    hit_table_task: ProteinHitTableTask
    gene_coords_task: Task
    gene_coords_reader: GeneCoordsReader
    fetcher: DenseWindowFetcher


@frozen(slots=True)
class ColocTaskGroup:
    coloc_task: TraitProteinColocTask
    display_task: ColocDisplayTableTask

    def terminal_tasks(self) -> list[Task]:
        return [self.display_task]


def generate_protein_coloc_tasks(
    base_name: str,
    source: ProteomicsColocSource,
    trait_task: Task,
    trait_model: TraitModel,
    config: ColocRunConfig,
    h4_threshold: float = 0.8,
) -> ColocTaskGroup:
    """
    base_name: prefix of the asset ids.
    trait_task: the trait GWAS, genome-reference harmonized to GRCh38.
    trait_model: quantitative or case-control, as coloc.abf needs.
    """
    prefix = f"{base_name}_{source.asset_infix}"
    coloc_task = TraitProteinColocTask.create(
        asset_id=f"{prefix}_coloc",
        trait_task=trait_task,
        trait_model=trait_model,
        hit_table_task=source.hit_table_task,
        gene_coords_task=source.gene_coords_task,
        gene_coords_reader=source.gene_coords_reader,
        fasta_task=UCSC_HG38_INDEXED_FASTA,
        panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
        panel_ancestry="nfe",
        fetcher=source.fetcher,
        config=config,
    )
    display_task = ColocDisplayTableTask.create(
        asset_id=f"{prefix}_coloc_display",
        coloc_task=coloc_task,
        h4_threshold=h4_threshold,
    )
    return ColocTaskGroup(coloc_task=coloc_task, display_task=display_task)


def generate_coloc_locus_plot_task(
    asset_id: str, group: ColocTaskGroup, pairs: tuple[ColocPair, ...]
) -> ColocLocusPlotTask:
    """Stacked locus plots of chosen tested pairs of one coloc run."""
    return ColocLocusPlotTask.create(
        asset_id=asset_id,
        coloc_task=group.coloc_task,
        gene_info_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED,
        pairs=pairs,
    )
