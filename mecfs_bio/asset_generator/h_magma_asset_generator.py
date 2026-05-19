"""
Asset generator for running H-MAGMA against GWAS summary statistics.

H-MAGMA (Sey et al. 2020) maps SNPs to genes using tissue-specific Hi-C
chromatin interaction data. The upstream project at
https://github.com/thewonlab/H-MAGMA/tree/master/Input_Files publishes six
pre-built ``.genes.annot`` files, one per tissue/cell type (adult brain, fetal
brain, cortical neurons, midbrain dopaminergic neurons, iPSC-derived
astrocytes, iPSC-derived neurons). This generator runs the MAGMA gene
analysis step against each of those six annotation files and produces a
gene-level Manhattan plot for each.

H-MAGMA's pre-built annotations replace the usual ``magma --annotate`` step,
so this generator skips :class:`MagmaAnnotateTask` and feeds the static
annotation directly into
:meth:`MagmaGeneAnalysisTask.create_with_prebuilt_annotation`.

The H-MAGMA annotation files are aligned to GRCh37/hg19 and key SNPs by RSID,
so the standard EUR build-37 1000 Genomes LD reference is used.
"""

from pathlib import PurePath

from attrs import frozen

from mecfs_bio.assets.executable.extracted.magma_binary_extracted import (
    MAGMA_1_1_BINARY_EXTRACTED,
)
from mecfs_bio.assets.reference_data.ensembl_biomart.gene_thesaurus import (
    GENE_THESAURUS,
)
from mecfs_bio.assets.reference_data.h_magma_annotations.raw.adult_brain_h_magma_annot import (
    ADULT_BRAIN_H_MAGMA_ANNOT_RAW,
)
from mecfs_bio.assets.reference_data.h_magma_annotations.raw.cortical_neuron_h_magma_annot import (
    CORTICAL_NEURON_H_MAGMA_ANNOT_RAW,
)
from mecfs_bio.assets.reference_data.h_magma_annotations.raw.fetal_brain_h_magma_annot import (
    FETAL_BRAIN_H_MAGMA_ANNOT_RAW,
)
from mecfs_bio.assets.reference_data.h_magma_annotations.raw.ipsc_derived_astro_h_magma_annot import (
    IPSC_DERIVED_ASTRO_H_MAGMA_ANNOT_RAW,
)
from mecfs_bio.assets.reference_data.h_magma_annotations.raw.ipsc_derived_neuro_h_magma_annot import (
    IPSC_DERIVED_NEURO_H_MAGMA_ANNOT_RAW,
)
from mecfs_bio.assets.reference_data.h_magma_annotations.raw.midbrain_da_h_magma_annot import (
    MIDBRAIN_DA_H_MAGMA_ANNOT_RAW,
)
from mecfs_bio.assets.reference_data.magma_ld_reference.magma_eur_build_37_1k_genomes_ref_extracted import (
    MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gene_manhattan_plot_task import (
    GeneManhattanPlotTask,
    MagmaGeneSource,
)
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild

# H-MAGMA annotations are built against GRCh37 / hg19.
_H_MAGMA_GENOME_BUILD: GenomeBuild = "19"

H_MAGMA_ANNOTATION_TASKS: list[tuple[str, Task]] = [
    ("adult_brain", ADULT_BRAIN_H_MAGMA_ANNOT_RAW),
    ("cortical_neuron", CORTICAL_NEURON_H_MAGMA_ANNOT_RAW),
    ("fetal_brain", FETAL_BRAIN_H_MAGMA_ANNOT_RAW),
    ("midbrain_da", MIDBRAIN_DA_H_MAGMA_ANNOT_RAW),
    ("ipsc_derived_astro", IPSC_DERIVED_ASTRO_H_MAGMA_ANNOT_RAW),
    ("ipsc_derived_neuro", IPSC_DERIVED_NEURO_H_MAGMA_ANNOT_RAW),
]


@frozen
class HMagmaTasksForAnnotation:
    """All tasks produced for a single H-MAGMA tissue annotation."""

    annotation_name: str
    gene_analysis_task: MagmaGeneAnalysisTask
    gene_manhattan_plot_task: GeneManhattanPlotTask


@frozen
class HMagmaTasks:
    """The aggregate result of running H-MAGMA against all six annotations."""

    p_value_task: MagmaSNPFileTask
    per_annotation: list[HMagmaTasksForAnnotation]

    def terminal_tasks(self) -> list[Task]:
        return [t.gene_manhattan_plot_task for t in self.per_annotation]


def generate_h_magma_tasks(
    base_name: str,
    gwas_parquet_with_rsids_task: Task,
    sample_size: int,
    pipes: list[DataProcessingPipe] | None = None,
) -> HMagmaTasks:
    """Generate one MAGMA gene analysis task and one gene-level Manhattan plot
    task per H-MAGMA tissue annotation (six in total).

    ``gwas_parquet_with_rsids_task`` must produce a parquet with the GWASLAB
    column names plus an RSID column (the standard input to MAGMA in this
    repository).
    """
    p_value_task = MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
        gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
        asset_id=base_name + "_h_magma_snp_p_values",
        pipes=pipes,
    )

    per_annotation: list[HMagmaTasksForAnnotation] = []
    for annotation_name, annotation_task in H_MAGMA_ANNOTATION_TASKS:
        gene_analysis_task = MagmaGeneAnalysisTask.create_with_prebuilt_annotation(
            asset_id=f"{base_name}_h_magma_{annotation_name}_gene_analysis",
            magma_annotation_task=annotation_task,
            magma_p_value_task=p_value_task,
            magma_binary_task=MAGMA_1_1_BINARY_EXTRACTED,
            magma_ld_ref_task=MAGMA_EUR_BUILD_37_1K_GENOMES_EXTRACTED,
            ld_ref_file_stem="g1000_eur",
            sample_size=sample_size,
            sub_dir_suffix=PurePath("h_magma") / annotation_name,
        )
        gene_manhattan_plot_task = GeneManhattanPlotTask.create(
            asset_id=f"{base_name}_h_magma_{annotation_name}_gene_manhattan_plot",
            source=MagmaGeneSource(
                magma_task=gene_analysis_task,
                gene_thesaurus_task=GENE_THESAURUS,
                genome_build=_H_MAGMA_GENOME_BUILD,
            ),
        )
        per_annotation.append(
            HMagmaTasksForAnnotation(
                annotation_name=annotation_name,
                gene_analysis_task=gene_analysis_task,
                gene_manhattan_plot_task=gene_manhattan_plot_task,
            )
        )

    return HMagmaTasks(
        p_value_task=p_value_task,
        per_annotation=per_annotation,
    )
