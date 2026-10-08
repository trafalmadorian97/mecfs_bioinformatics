"""
Where do the Python allocations go when building the UKB-PPP slim-protein Tasks?
Imports the index first (so its cost is excluded), then tracemallocs the generator module.
"""

import importlib
import tracemalloc

import mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index  # noqa: F401

tracemalloc.start(25)
snap0 = tracemalloc.take_snapshot()
mod = importlib.import_module(
    "mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_database_protein_files"
)
snap1 = tracemalloc.take_snapshot()
n = len(mod.HAPMAP_3_PPP_DATABASE.protein_tasks)
diff = snap1.compare_to(snap0, "lineno")
total = sum(s.size_diff for s in diff)
print(f"{n} tasks; traced python allocations {total / 2**20:.1f} MB ({total / n / 1024:.1f} KB/task)")
print("\n== by allocating line (mecfs_bio frames) ==")
for s in snap1.compare_to(snap0, "lineno")[:15]:
    print(f"{s.size_diff / 2**20:7.1f} MB  {s.traceback}")
print("\n== by traceback, innermost mecfs_bio frame ==")
from collections import defaultdict
agg = defaultdict(int)
for s in snap1.compare_to(snap0, "traceback"):
    frames = [f for f in s.traceback if "mecfs_bio" in f.filename]
    key = f"{frames[-1].filename.split('mecfs_bio/')[-1]}:{frames[-1].lineno}" if frames else "?"
    agg[key] += s.size_diff
for k, v in sorted(agg.items(), key=lambda kv: -kv[1])[:12]:
    print(f"{v / 2**20:7.1f} MB  {k}")
task = mod.HAPMAP_3_PPP_DATABASE.protein_tasks[0]
print("\nexample task:", type(task).__name__)
print(repr(task)[:1500])
