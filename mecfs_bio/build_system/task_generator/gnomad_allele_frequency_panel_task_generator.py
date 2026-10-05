"""
Builds a gnomAD release's allele-frequency panel together with its per-chromosome parts.
"""

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_allele_frequency_panel_task import (
    GnomadAlleleFrequencyPanelTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)


@frozen(slots=True)
class GnomadAlleleFrequencyPanelTasks:
    part_tasks: tuple[GnomadChromosomeAlleleFrequencyTask, ...]
    panel_task: GnomadAlleleFrequencyPanelTask

    def terminal_tasks(self) -> list[Task]:
        return [self.panel_task]


def generate_gnomad_allele_frequency_panel_tasks(
    asset_id: str, release: GnomadRelease, fasta_task: Task
) -> GnomadAlleleFrequencyPanelTasks:
    part_tasks = tuple(
        GnomadChromosomeAlleleFrequencyTask.create(
            release=release, chrom=chrom, fasta_task=fasta_task
        )
        for chrom in release.chromosomes
    )
    return GnomadAlleleFrequencyPanelTasks(
        part_tasks=part_tasks,
        panel_task=GnomadAlleleFrequencyPanelTask.create(
            asset_id=asset_id, release=release, part_tasks=part_tasks
        ),
    )
