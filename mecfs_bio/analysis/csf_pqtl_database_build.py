"""
Entry point that builds the full HapMap3 CSF pQTL database: the shared variant index
and all 7,008 per-aptamer slim beta/se/N files.

Each slim task downloads its aptamer's ~203 MB GWAS-SSF file into a scratch dir and
discards it after writing the ~6.8 MB aligned parquet, so the transient download
volume is large (~1.4 TB) but the durable output is ~47.5 GB. Route the output
subtree to a large disk via path_remap before running (see the CSF database plan).
"""

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_database_aptamer_files import (
    HAPMAP_3_CSF_DATABASE,
)


def go() -> None:
    """Build the shared index and every per-aptamer slim file for the HapMap3 CSF
    pQTL database."""
    DEFAULT_RUNNER.run(HAPMAP_3_CSF_DATABASE.terminal_tasks())


if __name__ == "__main__":
    go()
