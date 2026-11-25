"""
Is the PPP module's RSS cost polars' first-operation initialization rather than Tasks?
"""

import os

import polars as pl

from mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator import (
    EUR_DISCOVERY_PPP_MANIFEST_PATH,
)


def rss_mb() -> float:
    for line in open("/proc/self/status"):
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) / 1024
    raise AssertionError


print(f"POLARS_MAX_THREADS={os.environ.get('POLARS_MAX_THREADS')} thread_pool_size={pl.thread_pool_size()} cpus={os.cpu_count()}")
r0 = rss_mb()
pl.read_csv(EUR_DISCOVERY_PPP_MANIFEST_PATH)
r1 = rss_mb()
pl.read_csv(EUR_DISCOVERY_PPP_MANIFEST_PATH)
r2 = rss_mb()
import mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_database_protein_files  # noqa: E402,F401

r3 = rss_mb()
print(f"after imports: {r0:.0f} MB")
print(f"first pl.read_csv(manifest): +{r1 - r0:.1f} MB")
print(f"second pl.read_csv(manifest): +{r2 - r1:.1f} MB")
print(f"then import PPP module (index + 2940 tasks): +{r3 - r2:.1f} MB")
