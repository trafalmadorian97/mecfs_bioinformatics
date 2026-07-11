"""
Generates POPs gene-prioritization tasks using standard reference data.

Wires in the pinned POPs source and the shared munged feature collection so the
caller only supplies the GWAS-specific MAGMA gene-analysis task.
"""

from mecfs_bio.assets.reference_data.pops.features.pops_features_munged import (
    POPS_FEATURES_MUNGED,
)
from mecfs_bio.assets.reference_data.pops.source.pops_source_extracted import (
    POPS_SOURCE_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task_generator.pops_task_generator import PopsTaskGenerator


def concrete_pops_assets_generate(
    base_name: str,
    magma_gene_analysis_task: Task,
    pops_extra_args: tuple[str, ...] = (),
    low_mem: bool = False,
) -> PopsTaskGenerator:
    """
    Function to generate tasks that run POPs on top of a MAGMA gene analysis using
    standard reference data (the pinned POPs source and the shared munged POPs
    feature collection).

    magma_gene_analysis_task must be the ensembl MAGMA gene-analysis variant, since
    POPs matches genes on Ensembl gene IDs.

    Set low_mem=True to use the streaming kernel-ridge reimplementation, whose peak
    memory is independent of the selected-feature count. It is model-identical to
    stock POPs and lets high-feature traits run uncapped without OOMing the box.
    """
    return PopsTaskGenerator.create(
        base_name=base_name,
        pops_source_task=POPS_SOURCE_EXTRACTED,
        munged_features_task=POPS_FEATURES_MUNGED,
        magma_gene_analysis_task=magma_gene_analysis_task,
        pops_extra_args=pops_extra_args,
        low_mem=low_mem,
    )
