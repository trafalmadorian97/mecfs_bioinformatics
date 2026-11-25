"""
Readers for ldsc's chromosome-split LD-score files (<prefix><chr>.l2.ldscore.gz and
<prefix><chr>.l2.M_5_50), reproducing the parsing in ldsc's ldsc_parse module: chromosome
substitution, CHR/BP sort, first-occurrence SNP de-duplication, identical SNP columns across
files, and M_5_50 summed over the chromosomes present.

"""

from pathlib import Path
from typing import Sequence

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.sldsc.prepared_sldsc_reference import (
    BP_COL,
    CHR_COL,
    SNP_COL,
    LdScoreTable,
)

_MAX_CHROM = 22
_LD_SUFFIX = ".l2.ldscore.gz"
_M_SUFFIX = ".l2.M_5_50"
_KEY_COLS = (SNP_COL, CHR_COL, BP_COL)


def chromosome_file_prefix(prefix: str, chrom: int) -> str:
    """ldsc's sub_chr: append @ if absent, then substitute the chromosome number."""
    with_marker = prefix if "@" in prefix else prefix + "@"
    return with_marker.replace("@", str(chrom))


def present_chromosomes(root: Path, prefix: str) -> list[int]:
    return [
        c
        for c in range(1, _MAX_CHROM + 1)
        if (root / (chromosome_file_prefix(prefix, c) + _LD_SUFFIX)).is_file()
    ]


def read_ld_score_prefix(root: Path, prefix: str) -> pl.DataFrame:
    chroms = present_chromosomes(root, prefix)
    assert chroms, f"no LD-score files for prefix {prefix} under {root}"
    frames = [
        pl.read_csv(
            root / (chromosome_file_prefix(prefix, c) + _LD_SUFFIX), separator="\t"
        )
        for c in chroms
    ]
    frame = pl.concat(frames)
    value_cols = [c for c in frame.columns if c not in _KEY_COLS]
    return (
        frame.select(
            pl.col(SNP_COL).cast(pl.String),
            pl.col(CHR_COL).cast(pl.Int64),
            pl.col(BP_COL).cast(pl.Int64),
            *[pl.col(c).cast(pl.Float64) for c in value_cols],
        )
        .sort([CHR_COL, BP_COL], maintain_order=True)
        .unique(subset=SNP_COL, keep="first", maintain_order=True)
    )


def read_m_5_50_prefix(root: Path, prefix: str) -> np.ndarray:
    per_chrom = [
        [
            float(v)
            for v in (root / (chromosome_file_prefix(prefix, c) + _M_SUFFIX))
            .read_text()
            .splitlines()[0]
            .split()
        ]
        for c in present_chromosomes(root, prefix)
    ]
    return np.sum(np.array(per_chrom, dtype=np.float64), axis=0)


def read_ld_score_sources(
    root: Path, prefixes: Sequence[str], *, drop_zero_variance: bool
) -> LdScoreTable:
    """
    Read and side-by-side concatenate the LD scores of several file prefixes (ldsc's
    ldscore_fromlist). Column ids are ld_0, ld_1, ... in prefix-then-column order; sources
    record prefix:column. Streams one prefix at a time into a preallocated float32 matrix so
    that hundreds of prefixes do not all sit in memory as dataframes.

    drop_zero_variance: drop columns whose values have zero variance (ldsc's _check_variance,
    applied to the baseline only).
    """
    assert prefixes, "need at least one prefix"
    snps: pl.DataFrame | None = None
    blocks: list[np.ndarray] = []
    sources: list[str] = []
    m_values: list[np.ndarray] = []
    for prefix in prefixes:
        frame = read_ld_score_prefix(root, prefix)
        if snps is None:
            snps = frame.select(_KEY_COLS)
        assert frame[SNP_COL].equals(snps[SNP_COL]), (
            f"LD scores for {prefix} do not have the same SNPs as {prefixes[0]}"
        )
        value_cols = [c for c in frame.columns if c not in _KEY_COLS]
        m = read_m_5_50_prefix(root, prefix)
        assert m.shape[0] == len(value_cols), (
            f"{prefix}: {m.shape[0]} M_5_50 values for {len(value_cols)} LD columns"
        )
        values = frame.select(value_cols).to_numpy()
        keep = (
            np.var(values, axis=0, ddof=1) != 0
            if drop_zero_variance
            else np.ones(len(value_cols), dtype=bool)
        )
        blocks.append(values[:, keep].astype(np.float32))
        sources.extend(f"{prefix}:{c}" for c, k in zip(value_cols, keep) if k)
        m_values.append(m[keep])
    assert snps is not None
    matrix = np.concatenate(blocks, axis=1) if len(blocks) > 1 else blocks[0]
    return LdScoreTable(
        snps=snps,
        column_ids=tuple(f"ld_{i}" for i in range(matrix.shape[1])),
        sources=tuple(sources),
        matrix=matrix,
        m_5_50=np.concatenate(m_values),
    )
