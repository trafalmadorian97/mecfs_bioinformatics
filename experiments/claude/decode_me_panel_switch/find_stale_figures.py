"""
Find the DecodeME documentation figures whose content changed when DecodeME harmonization
switched to the gnomAD v2.1.1 nfe_nwe panel.

Scope is DecodeME's own pages: figures whose id starts with decode_me, excluding the
decode_me_minus_ GWAS-by-subtraction figures. The univariate MiXeR figures are also excluded,
since refitting MiXeR is slow; they keep their current versions. Other studies downstream of the panel are left
for a later update.

Every DecodeME figure task downstream of the gnomAD panel is exported (cached build, no forced rerun)
into a scratch figure directory. Each exported file is hashed and compared with the committed
figure manifest. A figure is stale if its hash differs from its manifest entry.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.find_stale_figures \
        2>&1 | tee experiments/claude/decode_me_panel_switch/find_stale_figures.log
"""

from pathlib import Path, PurePath

from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks
from mecfs_bio.figures.fig_constants import FIGURE_MANIFEST_PATH
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.generate_figures import FIGURE_EXPORTER
from mecfs_bio.figures.manifest import FigureManifest, scan_figure_dir

EXPORT_DIR = Path("experiments/claude/decode_me_panel_switch/figure_export")
PANEL_ID = GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED.asset_id
DECODE_ME_PREFIX = "decode_me"
SUBTRACTION_PREFIX = "decode_me_minus_"
UNIVARIATE_MIXER_PREFIX = "decode_me_univariate_mixer_"


def is_decode_me_figure(task: Task) -> bool:
    asset_id = str(task.asset_id)
    return asset_id.startswith(DECODE_ME_PREFIX) and not asset_id.startswith(
        (SUBTRACTION_PREFIX, UNIVARIATE_MIXER_PREFIX)
    )


def panel_dependent_decode_me_figure_tasks() -> list[Task]:
    return [
        task
        for task in ALL_FIGURE_TASKS
        if is_decode_me_figure(task) and PANEL_ID in find_tasks([task]).tasks
    ]


def main() -> None:
    tasks = panel_dependent_decode_me_figure_tasks()
    print(
        f"{len(tasks)} of {len(ALL_FIGURE_TASKS)} figure tasks are DecodeME and depend on "
        f"{PANEL_ID}",
        flush=True,
    )
    FIGURE_EXPORTER.export(tasks, fig_dir=EXPORT_DIR)
    manifest = FigureManifest.load(FIGURE_MANIFEST_PATH).figures
    exported = scan_figure_dir(EXPORT_DIR).figures
    changed: list[PurePath] = []
    unchanged: list[PurePath] = []
    for path, sha in sorted(exported.items()):
        if path not in manifest:
            print(f"NOT IN MANIFEST: {path}")
        elif manifest[path] == sha:
            unchanged.append(path)
        else:
            changed.append(path)
    print(f"\n{len(unchanged)} unchanged:")
    for path in unchanged:
        print(f"  {path}")
    print(f"\n{len(changed)} changed:")
    for path in changed:
        print(f"  {path}")


if __name__ == "__main__":
    main()
