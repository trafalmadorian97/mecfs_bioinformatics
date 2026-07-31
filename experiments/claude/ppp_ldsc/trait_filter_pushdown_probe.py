"""
Check that narrowing the trait to the context variants really does keep a large trait off the
peak memory of an rg run.

align_trait_to_context filters the trait to the LDSC context's rsIDs while it is still a lazy
frame, on the expectation that polars applies the predicate as it scans rather than after loading
everything. That expectation is the whole reason the change was made, so it is worth measuring
rather than assuming.

Peak resident memory is a high-water mark that never falls, so the variants cannot be compared
within one process: each is measured in a subprocess of its own. The variants are the eager read
the task used to do, the lazy filtered read with the default in-memory engine, and the same read
on the streaming engine.

Run: pixi r python experiments/claude/ppp_ldsc/trait_filter_pushdown_probe.py
"""

import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np
import polars as pl

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

# A trait much larger than the context, standing in for whole-genome-sequencing sumstats.
_TRAIT_ROW_COUNTS = (10_000_000, 20_000_000, 40_000_000)
# Roughly the size of the HapMap3 context the PPP tasks regress on.
_CONTEXT_ROWS = 1_000_000
_CHUNK_ROWS = 2_000_000

_VARIANTS = ("eager", "lazy_in_memory", "lazy_streaming")


def write_trait_parquet(path: Path, trait_rows: int) -> None:
    """Write the trait in chunks, so building the fixture does not itself dominate memory."""
    import pyarrow.parquet

    rng = np.random.default_rng(0)
    writer = None
    try:
        for start in range(0, trait_rows, _CHUNK_ROWS):
            size = min(_CHUNK_ROWS, trait_rows - start)
            table = pl.DataFrame(
                {
                    GWASLAB_RSID_COL: [f"rs{i}" for i in range(start, start + size)],
                    GWASLAB_EFFECT_ALLELE_COL: ["A"] * size,
                    GWASLAB_NON_EFFECT_ALLELE_COL: ["G"] * size,
                    GWASLAB_BETA_COL: rng.normal(0.0, 1.0, size),
                    GWASLAB_SE_COL: rng.uniform(0.01, 0.02, size),
                    GWASLAB_SAMPLE_SIZE_COLUMN: np.full(size, 50_000.0),
                }
            ).to_arrow()
            if writer is None:
                writer = pyarrow.parquet.ParquetWriter(path, table.schema)
            writer.write_table(table)
    finally:
        if writer is not None:
            writer.close()


def run_variant(variant: str, trait_path: Path) -> None:
    """Executed in a subprocess: read the trait one way and report this process's peak RSS."""
    import resource

    rsids = pl.Series(
        GWASLAB_RSID_COL, [f"rs{i}" for i in range(_CONTEXT_ROWS)]
    ).implode()
    lazy = pl.scan_parquet(trait_path)

    if variant == "eager":
        frame = pl.read_parquet(trait_path)
    elif variant == "lazy_in_memory":
        frame = lazy.filter(pl.col(GWASLAB_RSID_COL).is_in(rsids)).collect()
    elif variant == "lazy_streaming":
        frame = lazy.filter(pl.col(GWASLAB_RSID_COL).is_in(rsids)).collect(
            engine="streaming"
        )
    else:
        raise ValueError(f"unknown variant {variant}")

    peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    print(f"{variant}\t{frame.height}\t{peak_mb:.0f}")


def main() -> None:
    if len(sys.argv) == 3:
        run_variant(sys.argv[1], Path(sys.argv[2]))
        return

    print(f"context size: {_CONTEXT_ROWS:,} variants")
    print(
        f"\n{'trait rows':>12}{'file MB':>10}"
        + "".join(f"{variant + ' MB':>20}" for variant in _VARIANTS)
    )
    for trait_rows in _TRAIT_ROW_COUNTS:
        with TemporaryDirectory() as tmp:
            trait_path = Path(tmp) / "trait.parquet"
            write_trait_parquet(trait_path, trait_rows)
            peaks = []
            for variant in _VARIANTS:
                output = subprocess.run(
                    [sys.executable, __file__, variant, str(trait_path)],
                    capture_output=True,
                    text=True,
                    check=True,
                ).stdout.strip()
                _name, _rows, peak = output.split("\t")
                peaks.append(float(peak))
            print(
                f"{trait_rows:>12,}{trait_path.stat().st_size / 1e6:>10,.0f}"
                + "".join(f"{peak:>20,.0f}" for peak in peaks)
            )


if __name__ == "__main__":
    main()
