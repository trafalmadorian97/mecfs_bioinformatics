"""
Exclusive ("self") RSS cost of every module executed while importing all asset modules.

A meta-path hook wraps each loader's create_module/exec_module and measures RSS around
it, subtracting RSS growth attributable to nested imports.  Self costs are then summed
per top-level package (and per second-level package for mecfs_bio and gwaslab).
"""

import importlib
import importlib.abc
import sys
from collections import defaultdict
from pathlib import Path


def rss_mb() -> float:
    with open("/proc/self/status") as f:
        for line in f:
            if line.startswith("VmRSS:"):
                return int(line.split()[1]) / 1024
    raise AssertionError


self_cost: dict[str, float] = defaultdict(float)
stack: list[list[float]] = []  # per active module: [child_rss]


class Wrapped(importlib.abc.Loader):
    def __init__(self, inner, name):
        self.inner = inner
        self.name = name

    def _measure(self, fn, *args):
        start = rss_mb()
        stack.append([0.0])
        try:
            return fn(*args)
        finally:
            (child,) = stack.pop()
            total = rss_mb() - start
            self_cost[self.name] += total - child
            if stack:
                stack[-1][0] += total

    def create_module(self, spec):
        return self._measure(self.inner.create_module, spec)

    def exec_module(self, module):
        return self._measure(self.inner.exec_module, module)

    def __getattr__(self, item):
        return getattr(self.inner, item)


class Finder(importlib.abc.MetaPathFinder):
    def find_spec(self, fullname, path, target=None):
        for finder in sys.meta_path:
            if finder is self:
                continue
            spec = finder.find_spec(fullname, path, target) if hasattr(finder, "find_spec") else None
            if spec is not None:
                if spec.loader is not None and hasattr(spec.loader, "exec_module"):
                    spec.loader = Wrapped(spec.loader, fullname)
                return spec
        return None


def module_name(path: Path) -> str:
    parts = path.with_suffix("").parts
    if parts[-1] == "__init__":
        parts = parts[:-1]
    return ".".join(parts)


sys.meta_path.insert(0, Finder())
for path in sorted(Path("mecfs_bio/assets").rglob("*.py")):
    try:
        importlib.import_module(module_name(path))
    except Exception:
        pass
print(f"final RSS {rss_mb():.0f} MB; sum of self costs {sum(self_cost.values()):.0f} MB")


def group(name: str) -> str:
    parts = name.split(".")
    depth = 3 if parts[0] == "mecfs_bio" else 2 if parts[0] in {"gwaslab", "scipy"} else 1
    return ".".join(parts[:depth])


by_group: dict[str, float] = defaultdict(float)
for name, mb in self_cost.items():
    by_group[group(name)] += mb
print("\n== self RSS by package ==")
for g, mb in sorted(by_group.items(), key=lambda kv: -kv[1])[:30]:
    print(f"{mb:8.1f} MB  {g}")
print("\n== top individual modules (self) ==")
for name, mb in sorted(self_cost.items(), key=lambda kv: -kv[1])[:25]:
    print(f"{mb:8.1f} MB  {name}")
