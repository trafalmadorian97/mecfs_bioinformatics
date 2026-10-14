"""
Rebuild the DecodeME outputs recorded by record_decode_me_baseline, now harmonized against
the gnomAD v2.1.1 nfe_nwe panel.

The targets are exactly the assets the baseline recorded as built, so the comparison
(compare_decode_me_to_baseline) is like-for-like; tasks the baseline never built are not
built here either. The harmonization tasks' panel dependency changed, so the build system
rebuilds them and everything downstream.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.rebuild_decode_me \
        2>&1 | tee experiments/claude/decode_me_panel_switch/rebuild_decode_me.log
"""

import json

from experiments.claude.decode_me_panel_switch.record_decode_me_baseline import (
    BASELINE_DIR,
    MANIFEST_NAME,
    root_modules,
    root_tasks,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.tasks.simple_tasks import find_tasks


def baseline_built_ids() -> set[str]:
    records = json.loads((BASELINE_DIR / MANIFEST_NAME).read_text())
    return {r["asset_id"] for r in records if r["exists"]}


def main() -> None:
    built = baseline_built_ids()
    tasks = find_tasks(root_tasks(root_modules())).tasks
    targets: list[Task] = [
        task for asset_id, task in tasks.items() if str(asset_id) in built
    ]
    missing = built - {str(t.asset_id) for t in targets}
    assert not missing, f"baseline assets no longer in the task graphs: {sorted(missing)}"
    print(f"rebuilding {len(targets)} targets")
    DEFAULT_RUNNER.run(targets)
    print("rebuild finished")


if __name__ == "__main__":
    main()
