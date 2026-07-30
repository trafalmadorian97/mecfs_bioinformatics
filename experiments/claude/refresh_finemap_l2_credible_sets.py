"""
Force-rebuild and republish the DecodeME SUSIE fine-mapping figures affected by
the ``max_credible_sets`` fix in ``fine_mapping_asset_generator.py``.

Background
----------
The "L=2" fine-mapping task (``susie_finemap_2_credible_set_task``) was
previously created with ``max_credible_sets=1``, making it identical to the
"L=1" run.  That has been corrected to ``max_credible_sets=2``.

The verifying-trace rebuilder hashes *asset files*, not Python task
parameters, so changing ``max_credible_sets`` in the generator does not by
itself invalidate the cached SUSIE assets.  This script therefore:

1. Force-rebuilds the L=2 tasks via ``must_rebuild_transitive``.  The rebuilder
   additionally rebuilds every task that transitively depends on them: the L=2
   stackplots and the 4-run UpSet plots (which include the L=2 run as one of
   their inputs).
2. Re-exports the affected figures into ``docs/_figs``.  ``invoke
   publish-figures`` cannot be used to refresh *changed* figures --- its
   ``pull``/``generate_new_figures`` steps only fetch figures that are missing
   from disk, so it would just restore the stale copies from the release.  We
   therefore drive ``FIGURE_EXPORTER`` directly, which overwrites the figure
   directory with the freshly built assets.
3. Pushes: rehash ``docs/_figs``, update the committed manifest, and upload any
   new content-addressed blobs to the GitHub release.

After running this, commit ``mecfs_bio/figures/figures_manifest.json``.
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr1_174_128_548_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr6_97_505_620_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr15_54_925_638_locus_plalindindromes import (
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr17_50_237_377_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.fine_mapping.with_palindromes.susie_finemap_decode_me_37_chr20_47_653_230_locus_palindromes import (
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES,
)
from mecfs_bio.figures.fig_constants import FIGURE_DIRECTORY
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.generate_figures import FIGURE_EXPORTER
from mecfs_bio.figures.key_scripts.push_figures import push_figures

FINE_MAP_GROUPS = [
    DECODE_ME_GWAS_37_CHR1_174_128_548_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR6_97_505_620_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR_15_54_925_638_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR17_50_237_377_FINEMAP_PALINDROMES,
    DECODE_ME_GWAS_37_CHR20_47_653_000_FINEMAP_PALNDROMES,
]

# The fine-mapping figure tasks that appear in ALL_FIGURE_TASKS. We re-export
# exactly these (re-exporting the unchanged L=10 stackplots is harmless --- the
# push step only uploads blobs whose content hash actually changed).
FINE_MAP_FIGURE_TASKS = [
    task
    for group in FINE_MAP_GROUPS
    for task in (
        group.susie_stackplot_task,
        group.susie_finemap_2_credible_set_plot,
        group.upset_plot_task,
        group.upset_plot_task_pip001,
    )
    if task in ALL_FIGURE_TASKS
]


def refresh_finemap_l2_credible_sets():
    """Rebuild the L=2 SUSIE assets, then re-export and republish the figures."""
    # 1. Force-rebuild the L=2 tasks (and everything downstream of them).
    targets = [task for group in FINE_MAP_GROUPS for task in group.terminal_tasks()]
    must_rebuild_transitive = [
        group.susie_finemap_2_credible_set_task for group in FINE_MAP_GROUPS
    ]
    DEFAULT_RUNNER.run(
        targets,
        incremental_save=True,
        must_rebuild_transitive=must_rebuild_transitive,
    )

    # 2. Overwrite the stale figures in docs/_figs with the rebuilt assets.
    FIGURE_EXPORTER.export(FINE_MAP_FIGURE_TASKS, fig_dir=FIGURE_DIRECTORY)

    # 3. Update the manifest and upload any new blobs to the release.
    push_figures(figure_tasks=ALL_FIGURE_TASKS)


if __name__ == "__main__":
    refresh_finemap_l2_credible_sets()
