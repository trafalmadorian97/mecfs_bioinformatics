"""
Read-only check of whether already-built assets would verify under the verifying-trace
rebuilder, without running the runner (the runner has no dry run and would start rebuilding
anything that fails).

For every task in the transitive dependency closure of the requested tasks: resolve its path
with the default runner's meta_to_path, check the file or directory exists, recompute its
imohash trace, compare with the trace recorded in the info store, and compare each recorded
dependency trace with that dependency's current trace. The info store is never written.

Usage:
    pixi r python experiments/claude/csf_coloc/verify_traces.py ppp_coloc
    pixi r python experiments/claude/csf_coloc/verify_traces.py csf_database
"""

import sys
from collections.abc import Sequence
from typing import Literal

from attrs import define, field

from mecfs_bio.analysis.ppp_coloc_analysis import (
    DECODE_ME_PPP_COLOC,
    LIU_ET_AL_2023_IBD_PPP_COLOC,
    RA_SEROPOSITIVE_PPP_COLOC,
    WILLER_LDL_PPP_COLOC,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_database_aptamer_files import (
    HAPMAP_3_CSF_DATABASE,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.verifying_trace_info import (
    VerifyingTraceInfo,
)
from mecfs_bio.build_system.scheduler.asset_retrieval import get_asset_if_exists
from mecfs_bio.build_system.task.base_task import Task

Outcome = Literal[
    "verified", "missing_file", "no_trace", "value_mismatch", "dep_mismatch"
]
_OUTCOMES: tuple[Outcome, ...] = (
    "verified",
    "missing_file",
    "no_trace",
    "value_mismatch",
    "dep_mismatch",
)
_MAX_EXAMPLES = 20


@define(slots=True)
class VerifyReport:
    counts: dict[Outcome, int] = field(factory=lambda: {o: 0 for o in _OUTCOMES})
    examples: dict[Outcome, list[str]] = field(
        factory=lambda: {o: [] for o in _OUTCOMES}
    )

    def record(self, asset_id: AssetId, outcome: Outcome) -> None:
        self.counts[outcome] += 1
        if outcome != "verified" and len(self.examples[outcome]) < _MAX_EXAMPLES:
            self.examples[outcome].append(str(asset_id))

    def all_verified(self) -> bool:
        return all(self.counts[o] == 0 for o in _OUTCOMES if o != "verified")

    def render(self) -> str:
        lines = [f"{o}: {self.counts[o]}" for o in _OUTCOMES]
        for o in _OUTCOMES:
            if self.examples[o]:
                lines.append(f"examples {o}: {self.examples[o]}")
        return "\n".join(lines)


def _closure(tasks: Sequence[Task]) -> dict[AssetId, Task]:
    seen: dict[AssetId, Task] = {}
    stack = list(tasks)
    while stack:
        task = stack.pop()
        if task.asset_id in seen:
            continue
        seen[task.asset_id] = task
        stack.extend(task.deps)
    return seen


def verify_tasks(tasks: Sequence[Task]) -> VerifyReport:
    info = VerifyingTraceInfo.deserialize(DEFAULT_RUNNER.info_store)
    meta_to_path = DEFAULT_RUNNER.meta_to_path
    tracer = ImoHasher.with_xxhash_128()
    by_id = _closure(tasks)
    assets: dict[AssetId, Asset | None] = {}
    traces: dict[AssetId, str | None] = {}

    def asset_of(asset_id: AssetId) -> Asset | None:
        if asset_id not in assets:
            assets[asset_id] = get_asset_if_exists(by_id[asset_id].meta, meta_to_path)
        return assets[asset_id]

    def trace_of(asset_id: AssetId) -> str | None:
        if asset_id not in traces:
            asset = asset_of(asset_id) if asset_id in by_id else None
            traces[asset_id] = None if asset is None else tracer(asset)
        return traces[asset_id]

    report = VerifyReport()
    for asset_id in sorted(by_id):
        if asset_of(asset_id) is None:
            report.record(asset_id, "missing_file")
            continue
        if asset_id not in info.trace_store:
            report.record(asset_id, "no_trace")
            continue
        recorded_trace, recorded_deps = info.trace_store[asset_id]
        if trace_of(asset_id) != recorded_trace:
            report.record(asset_id, "value_mismatch")
            continue
        if any(trace_of(dep) != dep_trace for dep, dep_trace in recorded_deps):
            report.record(asset_id, "dep_mismatch")
            continue
        report.record(asset_id, "verified")
    return report


def _ppp_coloc_tasks() -> list[Task]:
    groups = [
        DECODE_ME_PPP_COLOC,
        LIU_ET_AL_2023_IBD_PPP_COLOC,
        RA_SEROPOSITIVE_PPP_COLOC,
        WILLER_LDL_PPP_COLOC,
    ]
    tasks: list[Task] = []
    for g in groups:
        tasks.extend([g.coloc_task, g.display_task, g.plot_task])
    tasks.append(groups[0].coloc_task.hit_table_task)
    for t in tasks:
        print(f"requested: {t.asset_id}")
    return tasks


def _csf_database_tasks() -> list[Task]:
    """The 7,008 slim aptamer files; their closure adds the index and the template aptamer."""
    return list(HAPMAP_3_CSF_DATABASE.aptamer_tasks)


def main() -> None:
    which = sys.argv[1]
    if which == "ppp_coloc":
        tasks = _ppp_coloc_tasks()
    elif which == "csf_database":
        tasks = _csf_database_tasks()
    else:
        raise ValueError(f"unknown mode {which}")
    report = verify_tasks(tasks)
    print(report.render())
    print("ALL VERIFIED" if report.all_verified() else "NOT ALL VERIFIED")


if __name__ == "__main__":
    main()
