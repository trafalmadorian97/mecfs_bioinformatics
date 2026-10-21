"""
Where does RSS go when load_sldsc_sumstats reads the DecodeME gwaslab pickle? Prints VmRSS after
each step of the load, with the whole frame converted to polars vs only the needed columns.

Usage:
    pixi r python experiments/claude/vectorized_sldsc/sumstats_load_rss.py <all|narrow>
"""

import gc
import sys
from pathlib import Path

import narwhals

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_sldsc import (
    DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS,
)
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePIfNeededPipe

_NEEDED = ["rsID", "EA", "NEA", "N", "BETA", "SE", "Z"]


def _rss_gb() -> float:
    line = next(
        x
        for x in Path("/proc/self/status").read_text().splitlines()
        if x.startswith("VmRSS:")
    )
    return int(line.split()[1]) / 1024**2


def main(mode: str) -> None:
    asset = DEFAULT_RUNNER.run([DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS])[
        DECODE_ME_GWAS_1_37_ANNOVAR_RSIDS_SUMSTATS.asset_id
    ]
    print(f"start {_rss_gb():.2f} GB")
    sumstats = read_sumstats(asset)
    data = sumstats.data
    print(
        f"pickle loaded {_rss_gb():.2f} GB; pandas frame {data.shape}, {data.memory_usage(deep=True).sum() / 1024**3:.2f} GB; columns {list(data.columns)}"
    )
    lazy = ComputePIfNeededPipe().process(narwhals.from_native(data).lazy())
    if mode == "narrow":
        lazy = lazy.select([c for c in _NEEDED if c in data.columns])
    frame = lazy.collect().to_polars()
    print(
        f"to_polars ({mode}) {_rss_gb():.2f} GB; polars frame {frame.estimated_size() / 1024**3:.2f} GB"
    )
    del sumstats, data, lazy
    gc.collect()
    print(f"pandas dropped {_rss_gb():.2f} GB")
    del frame
    gc.collect()
    print(f"polars dropped {_rss_gb():.2f} GB")


if __name__ == "__main__":
    main(sys.argv[1])
