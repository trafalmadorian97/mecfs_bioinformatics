"""
Break a cProfile dump of a full serial unit-test run into broad cost categories.

Each function's own time (tottime, excluding callees) is assigned to one category by
where the function lives (module file path) or, for C built-ins, by its name. Summing
own time over all functions gives the total profiled time, so the category shares add
up to 100%.

Categories of interest:
- filesystem: raw OS file calls (open/read/write/stat/mkdir/rename/unlink, shutil copies),
- hashing: md5/sha of file contents (the build system hashes every asset it produces),
- dataframe_io: parquet/csv/pickle readers and writers of polars, pandas and pyarrow
  (these mix disk access with encoding and decoding),
- r_rpy2: time inside R via rpy2,
- import: module importing (mostly test collection),
- typeguard: runtime type checking instrumentation,
- numeric, dataframe_compute, gwaslab, pytest, other.

Also prints the largest own-time functions per category, and cumulative time of a few
build-system entry points.

Usage: pixi r python experiments/claude/test_speed/profile_categories.py \
    experiments/claude/test_speed/suite.prof | tee experiments/claude/test_speed/profile_categories.log
"""

import pstats
import sys
from collections import defaultdict

FS_BUILTIN_MARKERS = (
    "posix.",
    "_io.open",
    "built-in method io.open",
    "of '_io.",
    "_io.FileIO",
    "_io.BufferedReader",
    "_io.BufferedWriter",
    "_io.TextIOWrapper",
    "_os.",
    "copyfile",
    "sendfile",
    "fsync",
)

PATH_CATEGORIES = [
    ("hashing", ("hashlib",)),
    ("r_rpy2", ("/rpy2/",)),
    ("typeguard", ("/typeguard/",)),
    ("import", ("<frozen importlib", "importlib/", "pkgutil", "zipimport")),
    ("filesystem", ("/shutil.py", "/pathlib", "/tempfile.py", "/os.py", "/genericpath.py", "/glob.py")),
    ("gwaslab", ("/gwaslab/",)),
    ("numeric", ("/numpy/", "/scipy/", "/statsmodels/", "/sklearn/", "/numba/")),
    ("dataframe_compute", ("/polars/", "/pandas/", "/pyarrow/", "/narwhals/")),
    ("pytest", ("/_pytest/", "/pluggy/", "/xdist/", "/testmon/")),
    ("matplotlib_plot", ("/matplotlib/", "/plotly/", "/seaborn/", "/PIL/")),
    ("mecfs_bio_python", ("/mecfs_bio/",)),
    ("tests_python", ("/test_mecfs_bio/",)),
]

DATAFRAME_IO_NAMES = (
    "read_parquet",
    "write_parquet",
    "scan_parquet",
    "sink_parquet",
    "read_csv",
    "write_csv",
    "scan_csv",
    "to_parquet",
    "to_csv",
    "read_table",
    "write_table",
    "read_ipc",
    "write_ipc",
    "collect",
    "pickle",
    "dump",
    "load",
)


def categorize(filename: str, name: str) -> str:
    is_builtin = filename == "~"
    if is_builtin:
        if "_hashlib" in name or "hashlib" in name:
            return "hashing"
        if "rpy2" in name or "_rinterface" in name or "cffi" in name:
            return "r_rpy2"
        if any(m in name for m in FS_BUILTIN_MARKERS):
            return "filesystem"
        if "marshal" in name or "_imp." in name:
            return "import"
        if "pickle" in name:
            return "dataframe_io"
        if "numpy" in name or "scipy" in name:
            return "numeric"
        if "polars" in name or "pyarrow" in name or "pandas" in name:
            return "dataframe_compute"
        if "sleep" in name or "select" in name or "poll" in name or "waitpid" in name:
            return "subprocess_wait"
        return "other_builtin"
    if ("/polars/" in filename or "/pandas/" in filename or "/pyarrow/" in filename) and any(
        n in name for n in DATAFRAME_IO_NAMES
    ):
        return "dataframe_io"
    if "/subprocess.py" in filename:
        return "subprocess_wait"
    for category, markers in PATH_CATEGORIES:
        if any(m in filename for m in markers):
            return category
    return "other_python"


def main(path: str) -> None:
    stats = pstats.Stats(path)
    raw = stats.stats  # type: ignore[attr-defined]
    totals: dict[str, float] = defaultdict(float)
    per_func: dict[str, list[tuple[float, str]]] = defaultdict(list)
    grand = 0.0
    for (filename, line, name), (_cc, _nc, tottime, _cumtime, _callers) in raw.items():
        category = categorize(filename, name)
        totals[category] += tottime
        per_func[category].append((tottime, f"{filename}:{line}({name})"))
        grand += tottime
    print(f"total own time: {grand:.1f}s")
    print(f"{'category':<20}{'seconds':>10}{'share':>8}")
    for category, secs in sorted(totals.items(), key=lambda kv: -kv[1]):
        print(f"{category:<20}{secs:>10.1f}{secs / grand:>8.1%}")
    for category in sorted(totals, key=lambda c: -totals[c]):
        print(f"\n--- top functions: {category} ---")
        for secs, label in sorted(per_func[category], reverse=True)[:8]:
            print(f"{secs:8.2f}s  {label[-140:]}")

    print("\n--- cumulative time of selected entry points ---")
    wanted = (
        "topological",
        "execute",
        "_rebuild",
        "hash_file",
        "md5",
        "materialize",
        "collect",
        "pytest_collection",
        "runtestprotocol",
    )
    rows = []
    for (filename, line, name), (_cc, nc, _tottime, cumtime, _callers) in raw.items():
        if name in wanted or any(name == w for w in wanted):
            rows.append((cumtime, nc, f"{filename}:{line}({name})"))
    for cumtime, nc, label in sorted(rows, reverse=True)[:25]:
        print(f"{cumtime:8.1f}s  calls={nc:<7} {label[-130:]}")


if __name__ == "__main__":
    main(sys.argv[1])
