"""
Attribute the resident memory of importing every mecfs_bio/assets module.

For each asset module, records the RSS growth caused by importing it and the top-level
packages that first appeared in sys.modules during that import.  The growth is charged
to those new packages, so heavy third-party libraries show up by name.
"""

import importlib
import resource
import sys
from collections import defaultdict
from pathlib import Path


def rss_mb() -> float:
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) / 1024
    raise AssertionError("no VmRSS")


def module_name(path: Path) -> str:
    parts = path.with_suffix("").parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


def top_levels() -> set[str]:
    return {name.split(".")[0] for name in sys.modules}


start = rss_mb()
print(f"baseline RSS: {start:.0f} MB")
charged: dict[str, float] = defaultdict(float)
per_module: list[tuple[float, str, list[str]]] = []
for path in sorted(Path("mecfs_bio/assets").rglob("*.py")):
    name = module_name(path)
    before_mods = top_levels()
    before = rss_mb()
    try:
        importlib.import_module(name)
    except Exception:
        pass
    delta = rss_mb() - before
    new = sorted(top_levels() - before_mods)
    per_module.append((delta, name, new))
    key = ", ".join(new) if new else "(mecfs_bio asset code / objects only)"
    charged[key] += delta

end = rss_mb()
print(f"final RSS: {end:.0f} MB (peak {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024:.0f} MB)")
print("\n== RSS growth by newly imported top-level packages ==")
for key, mb in sorted(charged.items(), key=lambda kv: -kv[1])[:25]:
    print(f"{mb:8.1f} MB  {key[:150]}")
print("\n== largest single-module deltas ==")
for delta, name, new in sorted(per_module, reverse=True)[:15]:
    print(f"{delta:8.1f} MB  {name}  new={new[:12]}")
