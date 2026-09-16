"""
V5: trust decision of every rsID-assignment chain's pre-harmonization table.

Imports every asset module, collects RSIDAssignmentTaskGroup instances, and for each group
whose upstream gwaslab Sumstats pickle is already materialized, builds the pre-harmonization
table and reports its trust evidence and decision. Groups whose inputs are not built are
listed as skipped, so nothing heavy is rebuilt. This tells us which sources take the
untrusted (stringent) path.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.survey_trust \
    2>&1 | tee experiments/claude/genome_reference_harmonization/survey_trust.log
"""

import importlib
import pkgutil

import mecfs_bio.assets
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.annovar_37_basic_rsid_assignment import (
    RSIDAssignmentTaskGroup,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import IndexedFasta
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import decide_trust
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

OPTIONS = GenomeReferenceHarmonizationOptions()


def _groups() -> list[RSIDAssignmentTaskGroup]:
    found: dict[str, RSIDAssignmentTaskGroup] = {}
    unimportable = 0
    for module_info in pkgutil.walk_packages(
        mecfs_bio.assets.__path__, "mecfs_bio.assets."
    ):
        try:
            module = importlib.import_module(module_info.name)
        except Exception as error:  # some asset modules fail at construction time
            unimportable += 1
            print(f"skip un-importable module {module_info.name}: {error!r}")
            continue
        for value in vars(module).values():
            if isinstance(value, RSIDAssignmentTaskGroup):
                found[value.harmonize_task.asset_id] = value
    print(f"found {len(found)} rsID-assignment groups; {unimportable} modules un-importable")
    return list(found.values())


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [UCSC_HG19_INDEXED_FASTA, THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES]
    )
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    panel_asset = assets[THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset) and isinstance(
        panel_asset, FileAsset
    )
    fasta = IndexedFasta.open(fasta_asset.path)
    for group in _groups():
        table_task = group.pre_harmonization_table_task
        assert isinstance(table_task, GwasLabSumstatsToTableTask)
        pickle_path = DEFAULT_RUNNER.meta_to_path(table_task.source_sumstats_task.meta)
        if not pickle_path.exists():
            print(f"SKIPPED (input not materialized): {group.harmonize_task.asset_id}")
            continue
        table_asset = DEFAULT_RUNNER.run([table_task])[table_task.asset_id]
        sumstats = scan_sumstats_as_polars(table_asset, table_task.meta, IdentityPipe())
        chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
        evidence = count_trust_evidence_genome_wide(
            sumstats, chromosomes, fasta, panel_asset.path, OPTIONS
        )
        trusted = decide_trust(evidence, OPTIONS)
        print(
            f"{group.harmonize_task.asset_id}: trusted={trusted} "
            f"eaf_present={evidence.eaf_present} {evidence.counts} "
            f"suspicious={evidence.suspicious}"
        )


if __name__ == "__main__":
    main()
