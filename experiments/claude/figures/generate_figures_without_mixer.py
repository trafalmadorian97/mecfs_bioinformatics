"""
Regenerate the documentation figures, skipping every figure that depends (transitively) on a
MiXeR task. MiXeR fits take hours, so this refreshes everything else without rerunning them;
the existing MiXeR figures in the figure directory are left as they are.

Usage:
    pixi r python experiments/claude/figures/generate_figures_without_mixer.py 2>&1 | tee experiments/claude/figures/generate_figures_without_mixer.log
"""

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.figures.figure_task_list import ALL_FIGURE_TASKS
from mecfs_bio.figures.key_scripts.generate_figures import generate_figures

_MIXER_MODULE_PREFIX = "mecfs_bio.build_system.task.mixer."


def _depends_on_mixer(task: Task, memo: dict[AssetId, bool]) -> bool:
    if task.asset_id not in memo:
        memo[task.asset_id] = type(task).__module__.startswith(
            _MIXER_MODULE_PREFIX
        ) or any(_depends_on_mixer(dep, memo) for dep in task.deps)
    return memo[task.asset_id]


def main() -> None:
    memo: dict[AssetId, bool] = {}
    kept = [t for t in ALL_FIGURE_TASKS if not _depends_on_mixer(t, memo)]
    skipped = [t for t in ALL_FIGURE_TASKS if _depends_on_mixer(t, memo)]
    print(f"skipping {len(skipped)} MiXeR figure(s):", flush=True)
    for t in skipped:
        print(f"  {t.asset_id}", flush=True)
    generate_figures(kept)


if __name__ == "__main__":
    main()
