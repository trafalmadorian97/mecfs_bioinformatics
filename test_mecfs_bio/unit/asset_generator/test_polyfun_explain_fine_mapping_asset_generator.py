from pathlib import PurePath

from mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator import (
    RUN_CONFIGS,
    PolyfunExplainOuterGroup,
    SharedFineMapInputs,
    build_explainability_groups,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask


def _shared() -> SharedFineMapInputs:
    harmonized = FakeTask(
        FilteredGWASDataMeta(
            id=AssetId("harmonized"),
            trait="mecfs",
            project="decodeme",
            sub_dir=PurePath("processed_gwas_data"),
        )
    )
    ld_labels = FakeTask(SimpleFileMeta("ld_labels"))
    ld_matrix = FakeTask(SimpleFileMeta("ld_matrix"))
    gene_info = FakeTask(SimpleFileMeta("genes"))
    return SharedFineMapInputs(
        base_name="mecfs_chr1_174",
        harmonized_sumstats_task=harmonized,
        ld_labels_task=ld_labels,
        ld_matrix_task=ld_matrix,
        gene_info_task=gene_info,
        effective_sample_size=10000,
        q_factor=100,
    )


def test_builds_four_groups_eight_susie_runs():
    groups = build_explainability_groups(_shared())
    assert len(groups) == len(RUN_CONFIGS) == 4
    susie_ids = []
    for g in groups:
        susie_ids += [g.susie_uniform.asset_id, g.susie_polyfun.asset_id]
    assert len(susie_ids) == 8
    assert len(set(susie_ids)) == 8  # all distinct


def test_polyfun_run_has_prior_uniform_does_not():
    group = build_explainability_groups(_shared())[0]
    assert group.susie_polyfun.prior_info is not None
    assert group.susie_uniform.prior_info is None


def test_each_group_wires_contrast_and_plot_to_its_own_pair():
    group = build_explainability_groups(_shared())[2]
    # The contrast/plot tasks explain THIS group's matched pair, so their deps
    # include the group's own uniform and polyfun runs.
    assert group.susie_uniform in group.contrast.deps
    assert group.susie_polyfun in group.contrast.deps
    assert group.contrast in group.plot.deps
    assert group.susie_uniform in group.plot.deps
    assert group.susie_polyfun in group.plot.deps


def test_outer_group_terminal_tasks_are_four_per_group():
    outer = PolyfunExplainOuterGroup(groups=build_explainability_groups(_shared()))
    # 4 configs * (uniform, polyfun, contrast, plot) = 16 terminal tasks.
    assert len(outer.terminal_tasks()) == 16
    assert len({t.asset_id for t in outer.terminal_tasks()}) == 16
