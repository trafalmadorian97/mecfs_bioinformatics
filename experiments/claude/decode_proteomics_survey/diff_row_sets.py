"""
Download one deCODE SMP file and compare its variant rows (Name column, in order) with a local
copy of another file of the same group:

  pixi r python experiments/claude/decode_proteomics_survey/diff_row_sets.py <seq_id> <reference.txt.gz> <scratch_dir>

Prints the rows each file has that the other lacks, and whether the shared rows keep one order.
"""

import io
import subprocess
import sys
from pathlib import Path

import polars as pl

from mecfs_bio.asset_generator.decode_somascan_manifest import load_decode_somascan_proteins
from mecfs_bio.build_system.task.decode_database.construct_decode_variant_index_task import (
    download_protein_file,
)
from mecfs_bio.build_system.wf.base_wf import make_wf


def _rows(path: Path) -> pl.DataFrame:
    text = subprocess.run(
        f"gzip -dc {path} | cut -f1,2,3,5,6,9", shell=True, check=True, capture_output=True
    ).stdout
    return pl.read_csv(
        io.BytesIO(text), separator="\t", infer_schema=False, quote_char=None
    ).with_row_index("row")


def main() -> None:
    seq_id, reference, scratch = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
    (protein,) = [p for p in load_decode_somascan_proteins("smp") if p.seq_id == seq_id]
    local = scratch / f"{seq_id}.txt.gz"
    if not local.exists():
        download_protein_file(protein, local, make_wf())
    a, b = _rows(local), _rows(reference)
    print(f"{seq_id}: {a.height} rows; reference: {b.height} rows")
    only_a = a.join(b, on="Name", how="anti")
    only_b = b.join(a, on="Name", how="anti")
    print(f"only in {seq_id}: {only_a.height}\n{only_a}")
    print(f"only in reference: {only_b.height}\n{only_b}")
    shared = a.join(b, on="Name", how="inner", suffix="_ref").sort("row")
    print(f"shared rows in the same order: {shared['row_ref'].is_sorted()}")
    print(f"duplicate Names: {seq_id} {a.height - a['Name'].n_unique()}, reference {b.height - b['Name'].n_unique()}")


if __name__ == "__main__":
    main()
