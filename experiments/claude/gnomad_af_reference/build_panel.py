"""
Validation step 2: build one allele-frequency panel and report its size and layout.

Usage (one panel per run):
    pixi r python -m experiments.claude.gnomad_af_reference.build_panel pan_ukbb \
        2>&1 | tee experiments/claude/gnomad_af_reference/build_panel_pan_ukbb.log
    pixi r python -m experiments.claude.gnomad_af_reference.build_panel gnomad_v2 \
        2>&1 | tee experiments/claude/gnomad_af_reference/build_panel_gnomad_v2.log
The gnomAD v2.1.1 build streams about 451 GiB and takes about 14 h; it resumes per
chromosome if interrupted. The Task logs (rows written, drops) are in the tee'd log.
"""

import sys

import polars as pl
import pyarrow.parquet as pq

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_allele_frequencies import (
    PAN_UKBB_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

PANELS: dict[str, Task] = {
    "pan_ukbb": PAN_UKBB_HG19_ALLELE_FREQUENCIES,
    "gnomad_v2": GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
}


def main(name: str) -> None:
    task = PANELS[name]
    asset = DEFAULT_RUNNER.run([task])[task.asset_id]
    assert isinstance(asset, FileAsset)
    metadata = pq.ParquetFile(asset.path).metadata
    print(f"{name}: {asset.path}")
    print(
        f"size {asset.path.stat().st_size / 2**30:.2f} GiB, rows {metadata.num_rows}, "
        f"row groups {metadata.num_row_groups}"
    )
    print(
        pl.scan_parquet(asset.path)
        .group_by(GWASLAB_CHROM_COL)
        .len()
        .sort(GWASLAB_CHROM_COL)
        .collect()
    )


if __name__ == "__main__":
    assert len(sys.argv) == 2 and sys.argv[1] in PANELS, (
        f"usage: build_panel.py {list(PANELS)}"
    )
    main(sys.argv[1])
