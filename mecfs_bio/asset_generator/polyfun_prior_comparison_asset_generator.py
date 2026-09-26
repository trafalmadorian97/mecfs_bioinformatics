"""Compare an external prior against an internal (trait-specific) prior at one
fine-mapping locus.

The external prior is fit outside the trait being fine-mapped (e.g. the
precomputed PolyFun prior); the internal prior is fit on the trait's own
sumstats (e.g. the L2-regularized S-LDSC prior). Both are run through
generate_assets_polyfun_explain_fine_map on the same locus, and the generators
here consume the two resulting PolyfunExplainOuterGroups:

- generate_polyfun_prior_comparison_group builds the comparison tasks for one
  run config from a uniform SUSIE run plus each prior's SUSIE run and contrast
  task.
- generate_assets_polyfun_prior_comparison pairs the two outer groups' run
  configs by label and builds one comparison group per config.
"""

from functools import cached_property
from pathlib import PurePath
from typing import Mapping

from attrs import frozen

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    PolyfunExplainGroup,
    PolyfunExplainOuterGroup,
)
from mecfs_bio.assets.reference_data.genetic_map.genetic_map_hg19 import (
    GENETIC_MAP_HG19,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.raw.magma_ensembl_gene_location_reference_data_build_37 import (
    MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    SecondaryPositionFromSnpid,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_plot_task import (
    PLOT_PNG_FILENAME,
    PLOT_SVG_FILENAME,
    PolyfunExplainPlotTask,
    PriorRun,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_table_task import (
    PolyfunPriorComparisonTableTask,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_comparison_variants import (
    PriorComparisonRuns,
)
from mecfs_bio.build_system.task.polyfun_explain.polyfun_prior_lift_scatter_task import (
    PolyfunPriorLiftScatterTask,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild

EXTERNAL_PRIOR_LABEL = "external prior"
INTERNAL_PRIOR_LABEL = "internal prior"


@frozen(slots=True)
class PolyfunPriorComparisonGroup:
    """The external-vs-internal prior comparison tasks for one run config."""

    label: str
    plot: Task
    # The plot's png and svg copied out of its directory as standalone FileAssets,
    # so docs can include either figure format on its own.
    plot_png: Task
    plot_svg: Task
    # Credible-set numbers, PIPs, and prior lifts of all three runs, one row per
    # credible-set variant of any run.
    table: Task
    # Interactive (html) scatter of the two priors' lifts over every locus
    # variant.
    lift_scatter: Task


@frozen(slots=True)
class PolyfunPriorComparisonOuterGroup:
    """One comparison group per run config for a single locus."""

    groups: list[PolyfunPriorComparisonGroup]

    def terminal_tasks(self) -> list[Task]:
        return [
            task
            for g in self.groups
            for task in (g.plot_png, g.plot_svg, g.table, g.lift_scatter)
        ]

    @cached_property
    def groups_by_label(self) -> Mapping[str, PolyfunPriorComparisonGroup]:
        return {group.label: group for group in self.groups}


def generate_polyfun_prior_comparison_group(
    base_name: str,
    label: str,
    susie_uniform: Task,
    external_prior_susie: Task,
    external_prior_contrast: Task,
    internal_prior_susie: Task,
    internal_prior_contrast: Task,
    gene_info_task: Task = MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    genome_build: GenomeBuild = "19",
    secondary_position: SecondaryPositionFromSnpid | None = None,
) -> PolyfunPriorComparisonGroup:
    """Build the comparison tasks for one run config: a stacked plot with PIP
    rows for the uniform, external-prior, and internal-prior runs, a table of
    the three runs' credible sets, PIPs, and prior lifts, and an interactive
    scatter of the two priors' lifts.

    secondary_position, when given, adds a build-labelled secondary position
    column (e.g. pos_hg38) to the table, parsed from the uniform run's SNPIDs."""
    stem = f"{base_name}_{label}"
    plot = PolyfunExplainPlotTask.create(
        asset_id=f"{stem}_prior_comparison_plot",
        susie_uniform_task=susie_uniform,
        prior_runs=(
            PriorRun(
                label=EXTERNAL_PRIOR_LABEL,
                susie_task=external_prior_susie,
                contrast_task=external_prior_contrast,
            ),
            PriorRun(
                label=INTERNAL_PRIOR_LABEL,
                susie_task=internal_prior_susie,
                contrast_task=internal_prior_contrast,
            ),
        ),
        gene_info_task=gene_info_task,
        genetic_map_task=GENETIC_MAP_HG19,
        genome_build=genome_build,
        gene_info_pipe=IdentityPipe(),
    )
    plot_png = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=f"{stem}_prior_comparison_plot_png",
        source_directory_task=plot,
        path_inside_directory=PurePath(PLOT_PNG_FILENAME),
        extension=".png",
    )
    plot_svg = CopyFileFromDirectoryTask.create_from_result_plot(
        asset_id=f"{stem}_prior_comparison_plot_svg",
        source_directory_task=plot,
        path_inside_directory=PurePath(PLOT_SVG_FILENAME),
        extension=".svg",
    )
    runs = PriorComparisonRuns(
        susie_uniform_task=susie_uniform,
        external_prior_susie_task=external_prior_susie,
        external_prior_contrast_task=external_prior_contrast,
        internal_prior_susie_task=internal_prior_susie,
        internal_prior_contrast_task=internal_prior_contrast,
    )
    table = PolyfunPriorComparisonTableTask.create(
        asset_id=f"{stem}_prior_comparison_table",
        runs=runs,
        secondary_position=secondary_position,
    )
    lift_scatter = PolyfunPriorLiftScatterTask.create(
        asset_id=f"{stem}_prior_lift_scatter",
        runs=runs,
        secondary_position=secondary_position,
    )
    return PolyfunPriorComparisonGroup(
        label=label,
        plot=plot,
        plot_png=plot_png,
        plot_svg=plot_svg,
        table=table,
        lift_scatter=lift_scatter,
    )


def generate_assets_polyfun_prior_comparison(
    base_name: str,
    external_prior: PolyfunExplainOuterGroup,
    internal_prior: PolyfunExplainOuterGroup,
    gene_info_task: Task = MAGMA_ENSEMBL_GENE_LOCATION_REFERENCE_DATA_BUILD_37_RAW,
    genome_build: GenomeBuild = "19",
    secondary_position: SecondaryPositionFromSnpid | None = None,
) -> PolyfunPriorComparisonOuterGroup:
    """Build one comparison group per run config shared by the two outer groups,
    which must come from the same locus and inputs, differing only in prior.

    Each outer group carries its own uniform SUSIE run; the two are computed from
    identical inputs, so the external group's is used as the common baseline."""
    external_labels = set(external_prior.groups_by_label)
    internal_labels = set(internal_prior.groups_by_label)
    assert external_labels == internal_labels, (
        f"Run configs differ between the external ({sorted(external_labels)}) "
        f"and internal ({sorted(internal_labels)}) prior outer groups"
    )
    groups = [
        _comparison_group_for_config(
            base_name=base_name,
            external=external,
            internal=internal_prior.groups_by_label[external.label],
            gene_info_task=gene_info_task,
            genome_build=genome_build,
            secondary_position=secondary_position,
        )
        for external in external_prior.groups
    ]
    return PolyfunPriorComparisonOuterGroup(groups=groups)


def _comparison_group_for_config(
    base_name: str,
    external: PolyfunExplainGroup,
    internal: PolyfunExplainGroup,
    gene_info_task: Task,
    genome_build: GenomeBuild,
    secondary_position: SecondaryPositionFromSnpid | None,
) -> PolyfunPriorComparisonGroup:
    assert external.label == internal.label
    return generate_polyfun_prior_comparison_group(
        base_name=base_name,
        label=external.label,
        susie_uniform=external.susie_uniform,
        external_prior_susie=external.susie_polyfun,
        external_prior_contrast=external.contrast,
        internal_prior_susie=internal.susie_polyfun,
        internal_prior_contrast=internal.contrast,
        gene_info_task=gene_info_task,
        genome_build=genome_build,
        secondary_position=secondary_position,
    )
