"""
Synthetic LDSC-SEG-format reference data and sumstats for S-LDSC unit tests.

Three chromosomes, 1,000 reference SNPs each. The baseline has a total-LD column, two
annotations, and an all-zero annotation (exercises ldsc's zero-variance drop). One shared
control plus three cell-type targets; the weight LD scores omit 30 SNPs on chromosome 2 (as the
real weights omit the MHC region); Cell_A carries a planted enrichment. Values are printed
with 3 decimals, like the real files. The sumstats omit some reference SNPs, add SNPs absent
from the reference (so gwaslab's smart_merge takes its pd.merge path, as in production), and
include a few chi-square outliers removed by the max(0.001 N, 80) filter.
"""

import gzip
from pathlib import Path

import gwaslab as gl
import numpy as np
import pandas as pd
from attrs import frozen

FIXTURE_BASELINE_PATTERN = "baseline/baseline.@"
FIXTURE_WEIGHTS_PATTERN = "weights/weights.@"
FIXTURE_LDCTS_FILENAME = "Set.ldcts"
FIXTURE_LABELS = ("Cell_A", "Cell_B", "Cell_C")
FIXTURE_N = 50_000

_CHROMS = (1, 2, 3)
_SNPS_PER_CHROM = 1_000
_BASELINE_COLS = ("baseL2", "AnnotAL2", "AnnotBL2", "ZeroL2")
_BASELINE_M = (5000.0, 1000.0, 500.0, 0.0)
_CONTROL_M = 2000.0
_TARGET_M = (300.0, 400.0, 500.0)
_WEIGHTLESS_CHROM = 2
_WEIGHTLESS_ROWS = slice(100, 130)


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
        rsids = [f"rs{chrom}{i:05d}" for i in range(_SNPS_PER_CHROM)]
        bps = 10_000 + 100 * np.arange(_SNPS_PER_CHROM)
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
        control = rng.uniform(0, 20, _SNPS_PER_CHROM)
        _write_ld_file(
            path=cts_dir / f"cts/Set.control.{chrom}.l2.ldscore.gz",
            chrom=chrom,
            rsids=rsids,
            bps=bps,
            cols={"L2": control},
        )
        (cts_dir / f"cts/Set.control.{chrom}.l2.M_5_50").write_text(
            f"{_CONTROL_M / len(_CHROMS)}\n"
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
            {"rsID": r, "CHR": chrom, "POS": int(bp), "Z": float(zz)}
            for r, bp, zz in zip(rsids, bps, z)
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
        }
    )
    frame = pd.concat([frame, extra], ignore_index=True)
    frame.loc[[100, 200, 300], "Z"] = 12.0  # chi2 = 144 > max(0.001 * N, 80)
    se = np.full(len(frame), 0.01)
    frame = frame.assign(
        BETA=frame["Z"] * se, SE=se, EA="A", NEA="G", P=0.5, N=FIXTURE_N
    ).drop(columns="Z")
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
