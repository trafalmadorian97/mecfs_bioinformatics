"""
Validation step 1: the panel-ancestry refactor leaves 1000 Genomes harmonizations unchanged.

Reads DecodeME's current harmonized table (built by the pre-refactor code), force-rebuilds it
with the current code, and compares the two as frames and as bytes. The build cache does not
track code, hence the forced rebuild.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.byte_identity_1000g \
        2>&1 | tee experiments/claude/gnomad_af_reference/byte_identity_1000g.log
"""

import hashlib
from pathlib import Path

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task


def _built_path(task: Task, rebuild: bool) -> Path:
    assets = DEFAULT_RUNNER.run(
        [task], must_rebuild_transitive=[task] if rebuild else ()
    )
    asset = assets[task.asset_id]
    assert isinstance(asset, FileAsset)
    return asset.path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    task = DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task
    before_path = _built_path(task, rebuild=False)
    before_frame = pl.read_parquet(before_path)
    before_hash = _sha256(before_path)
    after_path = _built_path(task, rebuild=True)
    print(f"rows before {before_frame.height}")
    print("frames equal:", before_frame.equals(pl.read_parquet(after_path)))
    print("bytes equal:", before_hash == _sha256(after_path))


if __name__ == "__main__":
    main()
