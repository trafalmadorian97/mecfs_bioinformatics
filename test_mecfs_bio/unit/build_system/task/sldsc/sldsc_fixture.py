"""
Synthetic LDSC-SEG-format reference data and sumstats for S-LDSC unit tests.

Three chromosomes, 1,000 reference SNPs each. The baseline has a total-LD column, two
annotations, and an all-zero annotation (exercises ldsc's zero-variance drop). One shared
control (two LD columns, like the 12-column Roadmap control of Corces_ATAC) plus three
cell-type targets; the weight LD scores omit 30 SNPs on chromosome 2 (as the
real weights omit the MHC region); Cell_A carries a planted enrichment. Values are printed
with 3 decimals, like the real files. The sumstats omit some reference SNPs, add SNPs absent
from the reference (so gwaslab's smart_merge takes its pd.merge path, as in production), and
include a few chi-square outliers removed by the max(0.001 N, 80) filter.

Reference SNPs are real HapMap3 SNPs (rsID, position, alleles from gwaslab's bundled build-19
HapMap3 list), because gwaslab's estimate_h2_cts_by_ldsc keeps only sumstats rows whose rsID is
in that list with matching alleles. Every third SNP reports its alleles swapped (kept), and every
25th SNP from the eighth reports an allele absent from HapMap3 (dropped).
"""

import gzip
from pathlib import Path

import gwaslab as gl
import numpy as np
import pandas as pd
import polars as pl
from attrs import frozen
from gwaslab.bd.bd_get_hapmap3 import _get_hapmap_full_polars

FIXTURE_BASELINE_PATTERN = "baseline/baseline.@"
FIXTURE_WEIGHTS_PATTERN = "weights/weights.@"
FIXTURE_LDCTS_FILENAME = "Set.ldcts"
FIXTURE_LABELS = ("Cell_A", "Cell_B", "Cell_C")
FIXTURE_N = 50_000

_CHROMS = (1, 2, 3)
_SNPS_PER_CHROM = 1_000
_BASELINE_COLS = ("baseL2", "AnnotAL2", "AnnotBL2", "ZeroL2")
_BASELINE_M = (5000.0, 1000.0, 500.0, 0.0)
_CONTROL_M = (2000.0, 1500.0)
_TARGET_M = (300.0, 400.0, 500.0)
_WEIGHTLESS_CHROM = 2
_WEIGHTLESS_ROWS = slice(100, 130)
_HAPMAP3_BUILD = "19"
_SWAPPED_ALLELE_PERIOD = 3
_MISMATCHED_ALLELE_PERIOD = 25
_MISMATCHED_ALLELE_OFFSET = 7
_NUCLEOTIDES = ("A", "C", "G", "T")


def _hapmap3_snps(chrom: int) -> pl.DataFrame:
    """The first _SNPS_PER_CHROM HapMap3 SNPs of a chromosome with a unique position."""
    return (
        _get_hapmap_full_polars(_HAPMAP3_BUILD, include_alleles=True)
        .filter(pl.col("CHR") == chrom)
        .filter(pl.col("POS").is_unique())
        .sort("POS")
        .head(_SNPS_PER_CHROM)
    )


def _reported_alleles(index: int, a1: str, a2: str) -> tuple[str, str]:
    if index % _MISMATCHED_ALLELE_PERIOD == _MISMATCHED_ALLELE_OFFSET:
        return next(n for n in _NUCLEOTIDES if n not in (a1, a2)), a2
    if index % _SWAPPED_ALLELE_PERIOD == 0:
        return a2, a1
    return a1, a2


@frozen(slots=True)
class SLDSCFixture:
    baseline_dir: Path
    weights_dir: Path
    cts_dir: Path
    sumstats_frame: pd.DataFrame


def _write_ld_file(
    path: Path,
    chrom: int,
    rsids: list[str],
    bps: np.ndarray,
    cols: dict[str, np.ndarray],
) -> None:
    with gzip.open(path, "wt") as f:
        f.write("\t".join(["CHR", "SNP", "BP", *cols]) + "\n")
        for i, rsid in enumerate(rsids):
            values = "\t".join(f"{cols[c][i]:.3f}" for c in cols)
            f.write(f"{chrom}\t{rsid}\t{bps[i]}\t{values}\n")


def write_sldsc_fixture(root: Path, *, seed: int = 0) -> SLDSCFixture:
    rng = np.random.default_rng(seed)
    baseline_dir = root / "baseline_root"
    weights_dir = root / "weights_root"
    cts_dir = root / "cts_root"
    for d in (baseline_dir / "baseline", weights_dir / "weights", cts_dir / "cts"):
        d.mkdir(parents=True)
    rows = []
    for chrom in _CHROMS:
        hapmap3 = _hapmap3_snps(chrom)
        rsids = hapmap3["rsid"].to_list()
        bps = hapmap3["POS"].to_numpy()
        alleles = [
            _reported_alleles(i, a1, a2)
            for i, (a1, a2) in enumerate(zip(hapmap3["A1"], hapmap3["A2"]))
        ]
        base = {
            "baseL2": rng.uniform(1, 50, _SNPS_PER_CHROM),
            "AnnotAL2": rng.uniform(0, 10, _SNPS_PER_CHROM),
            "AnnotBL2": rng.uniform(0, 5, _SNPS_PER_CHROM),
            "ZeroL2": np.zeros(_SNPS_PER_CHROM),
        }
        _write_ld_file(
            path=baseline_dir / f"baseline/baseline.{chrom}.l2.ldscore.gz",
            chrom=chrom,
            rsids=rsids,
            bps=bps,
            cols=base,
        )
        (baseline_dir / f"baseline/baseline.{chrom}.l2.M_5_50").write_text(
            "\t".join(str(m / len(_CHROMS)) for m in _BASELINE_M) + "\n"
        )
        weight_ld = rng.uniform(1, 40, _SNPS_PER_CHROM)
        has_weight = np.ones(_SNPS_PER_CHROM, dtype=bool)
        if chrom == _WEIGHTLESS_CHROM:
            has_weight[_WEIGHTLESS_ROWS] = False
        _write_ld_file(
            path=weights_dir / f"weights/weights.{chrom}.l2.ldscore.gz",
            chrom=chrom,
            rsids=[r for r, keep in zip(rsids, has_weight) if keep],
            bps=bps[has_weight],
            cols={"L2": weight_ld[has_weight]},
        )
        control = {
            "DNase.union": rng.uniform(0, 20, _SNPS_PER_CHROM),
            "DNase.sum": rng.uniform(0, 60, _SNPS_PER_CHROM),
        }
        _write_ld_file(
            path=cts_dir / f"cts/Set.control.{chrom}.l2.ldscore.gz",
            chrom=chrom,
            rsids=rsids,
            bps=bps,
            cols=control,
        )
        (cts_dir / f"cts/Set.control.{chrom}.l2.M_5_50").write_text(
            "\t".join(str(m / len(_CHROMS)) for m in _CONTROL_M) + "\n"
        )
        targets = [rng.uniform(0, 10, _SNPS_PER_CHROM) for _ in FIXTURE_LABELS]
        for j, target in enumerate(targets, start=1):
            _write_ld_file(
                path=cts_dir / f"cts/Set.{j}.{chrom}.l2.ldscore.gz",
                chrom=chrom,
                rsids=rsids,
                bps=bps,
                cols={"L2": target},
            )
            (cts_dir / f"cts/Set.{j}.{chrom}.l2.M_5_50").write_text(
                f"{_TARGET_M[j - 1] / len(_CHROMS)}\n"
            )
        expected_chi2 = 1.0 + 0.02 * base["baseL2"] + 0.08 * targets[0]
        z = rng.normal(0.0, 1.0, _SNPS_PER_CHROM) * np.sqrt(expected_chi2)
        rows.extend(
            {
                "rsID": r,
                "CHR": chrom,
                "POS": int(bp),
                "Z": float(zz),
                "EA": ea,
                "NEA": nea,
            }
            for r, bp, zz, (ea, nea) in zip(rsids, bps, z, alleles)
        )
    (cts_dir / FIXTURE_LDCTS_FILENAME).write_text(
        "\n".join(
            f"{label}\tcts/Set.{j}.,cts/Set.control."
            for j, label in enumerate(FIXTURE_LABELS, start=1)
        )
        + "\n"
    )
    frame = pd.DataFrame(rows)
    frame = frame.iloc[50:].reset_index(
        drop=True
    )  # some reference SNPs absent from the sumstats
    extra = pd.DataFrame(
        {
            "rsID": [f"rs_extra{i}" for i in range(20)],
            "CHR": 1,
            "POS": 9_000_000 + np.arange(20),
            "Z": rng.normal(0, 1, 20),
            "EA": "A",
            "NEA": "G",
        }
    )
    frame = pd.concat([frame, extra], ignore_index=True)
    frame.loc[[100, 200, 300], "Z"] = 12.0  # chi2 = 144 > max(0.001 * N, 80)
    se = np.full(len(frame), 0.01)
    frame = frame.assign(BETA=frame["Z"] * se, SE=se, P=0.5, N=FIXTURE_N).drop(
        columns="Z"
    )
    return SLDSCFixture(
        baseline_dir=baseline_dir,
        weights_dir=weights_dir,
        cts_dir=cts_dir,
        sumstats_frame=frame,
    )


def sumstats_from_frame(frame: pd.DataFrame) -> gl.Sumstats:
    return gl.Sumstats(
        frame,
        rsid="rsID",
        chrom="CHR",
        pos="POS",
        ea="EA",
        nea="NEA",
        beta="BETA",
        se="SE",
        p="P",
        n="N",
        build="19",
        verbose=False,
    )
