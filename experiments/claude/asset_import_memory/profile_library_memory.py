"""
Break down the RSS of the third-party stack loaded by the first asset import.

Prints, for each top-level library loaded by that import, (a) its standalone RSS cost in
a fresh interpreter and (b) its incremental cost when imported cumulatively in order of
standalone cost (largest first), so shared dependencies are charged once.
"""

import importlib
import subprocess
import sys

FIRST_ASSET = "mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.analysis.bellenguez_pp_rg"

RSS_SNIPPET = """
def rss():
    for line in open('/proc/self/status'):
        if line.startswith('VmRSS:'):
            return int(line.split()[1]) / 1024
"""


def run(code: str) -> str:
    return subprocess.run(
        [sys.executable, "-c", RSS_SNIPPET + code],
        capture_output=True, text=True, check=True,
    ).stdout.strip().splitlines()[-1]


listed = run(
    f"import sys; b=set(sys.modules); import {FIRST_ASSET}; "
    "print(','.join(sorted({m.split('.')[0] for m in sys.modules} - {m.split('.')[0] for m in b})))"
)
candidates = [
    m for m in listed.split(",")
    if not m.startswith("_") and m not in sys.stdlib_module_names and m != "mecfs_bio"
]
print(f"baseline interpreter RSS: {run('print(rss())')} MB")
print(f"asset import RSS: {run(f'import {FIRST_ASSET}; print(rss())')} MB")

standalone = {}
for m in candidates:
    try:
        standalone[m] = float(run(f"b=rss(); import {m}; print(rss()-b)"))
    except subprocess.CalledProcessError:
        pass
order = sorted(standalone, key=lambda m: -standalone[m])

cumulative = run(
    "import importlib, json; out={}\n"
    f"for m in {order!r}:\n"
    "    b=rss(); importlib.import_module(m); out[m]=rss()-b\n"
    "print(json.dumps(out))"
)
import json
incremental = json.loads(cumulative)
print(f"\n{'library':28s} {'standalone MB':>14s} {'incremental MB':>15s}")
for m in order:
    if standalone[m] >= 3 or incremental[m] >= 3:
        print(f"{m:28s} {standalone[m]:14.1f} {incremental[m]:15.1f}")
print(f"{'TOTAL incremental':28s} {'':14s} {sum(incremental.values()):15.1f}")
