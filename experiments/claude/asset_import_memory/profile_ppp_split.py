"""
Split the PPP module's RSS between its imports and its Task generation, with polars warm.
"""

import gc
import importlib

import polars as pl

import mecfs_bio.asset_generator.ukbb_ppp_slim_protein_asset_generator as gen
import mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index as idx


def rss_mb() -> float:
    for line in open("/proc/self/status"):
        if line.startswith("VmRSS:"):
            return int(line.split()[1]) / 1024
    raise AssertionError


pl.read_csv(gen.EUR_DISCOVERY_PPP_MANIFEST_PATH)
r0 = rss_mb()
manifest = pl.read_csv(gen.EUR_DISCOVERY_PPP_MANIFEST_PATH)
rows = list(manifest.iter_rows(named=True))
r1 = rss_mb()
protein = gen._protein_file_from_row(rows[0])
task = gen.BuildSlimProteinParquetTask.create(
    index_task=idx.HAPMAP_3_PPP_DATABASE_INDEX, protein=protein, asset_id="x",
    index_name="hapmap_3", include_sample_size=False,
)
r2 = rss_mb()
tasks = [
    gen.BuildSlimProteinParquetTask.create(
        index_task=idx.HAPMAP_3_PPP_DATABASE_INDEX, protein=gen._protein_file_from_row(r),
        asset_id=f"x{i}", index_name="hapmap_3", include_sample_size=False,
    )
    for i, r in enumerate(rows[:500])
]
r3 = rss_mb()
print(f"warm baseline {r0:.0f} MB")
print(f"manifest rows: +{r1 - r0:.1f} MB")
print(f"first task.create: +{r2 - r1:.1f} MB")
print(f"next 500 tasks: +{r3 - r2:.1f} MB ({(r3 - r2) / 500 * 1024:.1f} KB/task)")
