# PolyFun Approach 2 (L2 S-LDSC trait-specific prior) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Estimate a per-SNP heritability (`snpvar`) prior *from the GWAS being
fine-mapped* via an L2-regularized S-LDSC over the baseline-LF annotations, and
expose it as a drop-in alternative to the precomputed Approach-1 prior in the
SUSIE fine-mapping generator.

**Architecture:** A new estimator Task streams the baseline-LF annotation
LD-scores by chromosome, accumulates weighted ridge cross-product sufficient
statistics (never materializing the ~19.5M x 187 design matrix), fits annotation
coefficients `tau` nested within an odd/even chromosome split, and scores a dense
per-variant `snpvar`. The ridge linear algebra is a pure-numpy submodule shared
with the existing `RidgeAnnotationWeightsTask`. The prior plugs into the existing
`PriorInfo` seam.

**Tech Stack:** Python, polars (new dataframe code), numpy, pyarrow, attrs
(`@frozen`), the repo build-system Task framework, pixi (`pixi r ...`).

**Spec:** `experiments/claude/design_specs/2026-09-20-polyfun-approach-2-l2-sldsc-prior-design.md`
(read it — the plan argues from it). Schema spike:
`experiments/claude/polyfun_approach2_ldscore_schema_spike/`.

## Global Constraints

- **Join key:** every per-variant join uses the exact four-part key
  `(chrom, pos, nea, ea)` — gwaslab columns `CHR, POS, NEA, EA`
  (`GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_NON_EFFECT_ALLELE_COL,
  GWASLAB_EFFECT_ALLELE_COL`). No rsid, unordered-allele, or position-only join.
  baseline-LF members are REF-oriented: their `A1 == REF == nea`,
  `A2 == ALT == ea`. If any input lacks the four-part key, STOP.
- **PolyFun column constants:** every read, rename, select, or fixture of a
  PolyFun-provided dataframe (annotation members, LD-score members, the annotation
  matrix, the precomputed prior) names its key columns through
  `mecfs_bio/constants/polyfun_constants.py` (created in Task 1 Step 0):
  `POLYFUN_CHR_COL, POLYFUN_BP_COL, POLYFUN_SNP_COL, POLYFUN_A1_COL,
  POLYFUN_A2_COL`, plus `POLYFUN_TO_GWASLAB_KEY_RENAME`, the single mapping
  used to turn a PolyFun key into the gwaslab four-part key. No `"CHR"`/`"BP"`/
  `"A1"`/`"A2"` literals in new code or tests.
- **Memory budget:** 16GB. Never materialize the dense design matrix; stream one
  chromosome at a time. Large extracted assets are `path_remap` candidates.
- **No leakage (spec decision 8):** nothing that scores a chromosome may be fit
  on data from its own parity — weights (`h2bar`), ridge lambda, and coefficients
  are all fit per-parity on the opposite half.
- **Conventions:** polars over pandas; `Path`/`PurePath` in memory; column-name
  constants (no repeated literals); named kwargs for same-typed params; inject
  dependencies (no monkeypatch); `@frozen` attrs with `__attrs_post_init__`
  shape/dtype asserts on numpy fields; imports at top of file; `.create()` derives
  `meta` from a dependency; helper free functions over methods; run
  `pixi r invoke green` (capture to a logfile) after each task.

---

### Task 1: baseline-LF annotation LD-score members (stream-extract + asset)

**Files:**
- Create: `mecfs_bio/constants/polyfun_constants.py`
- Modify (constants only): `mecfs_bio/build_system/task/annotation_weights/build_baseline_lf_annotation_parquet_task.py`,
  `mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py`,
  `mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py`
- Modify: `mecfs_bio/build_system/task/annotation_weights/stream_extract_annotation_parquets_task.py`
- Create: `mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_ldscores.py`
- Test: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_stream_extract_annotation_parquets_task.py` (extend)

**Interfaces:**
- Produces: `BASELINE_LF_ANNOTATION_LDSCORE_MEMBERS: Task` — a `DirectoryAsset` of
  `baselineLF2.2.UKB.<chr>.l2.ldscore.parquet` members, each with columns
  `CHR, SNP, BP, A1, A2` + 187 annotation-named LD-score columns.

The existing `StreamExtractAnnotationParquetsTask` hardcodes the annot-member
regex and destination name. Generalize it with two injected, defaulted fields so a
second instantiation extracts the LD-score members, then add the asset.

- [ ] **Step 0: Add the PolyFun column constants and adopt them**

Create `mecfs_bio/constants/polyfun_constants.py`:

```python
"""Column names used by PolyFun-provided dataframes.

The baseline-LF annotation members, LD-score members, derived annotation matrix,
and precomputed prior all key variants by CHR, BP, SNP, A1, A2. A1 and A2 are
REF-oriented: A1 is the hg19 reference allele (the gwaslab non-effect allele) and
A2 is the alternate allele (the gwaslab effect allele).
"""

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

POLYFUN_CHR_COL = "CHR"
POLYFUN_BP_COL = "BP"
POLYFUN_SNP_COL = "SNP"
POLYFUN_A1_COL = "A1"
POLYFUN_A2_COL = "A2"

# Renames a PolyFun variant key to the gwaslab four-part join key.
POLYFUN_TO_GWASLAB_KEY_RENAME: dict[str, str] = {
    POLYFUN_CHR_COL: GWASLAB_CHROM_COL,
    POLYFUN_BP_COL: GWASLAB_POS_COL,
    POLYFUN_A1_COL: GWASLAB_NON_EFFECT_ALLELE_COL,
    POLYFUN_A2_COL: GWASLAB_EFFECT_ALLELE_COL,
}
```

Then replace the duplicated private copies with these shared constants (a pure
rename, no behavior change):
- `build_baseline_lf_annotation_parquet_task.py`: `_CHR_COL/_BP_COL/_A1_COL/_A2_COL`
  and the literals in `ANNOT_KEY_COLUMNS`.
- `ridge_annotation_weights_task.py`: `_CHR_COL/_BP_COL/_A1_COL/_A2_COL`.
- `susie_r_finemap_task.py`: the `PriorInfo` defaults `prior_chr_col`,
  `prior_bp_cp`, `prior_a1_col`, `prior_a2_col`.
- the existing `_annot_frame` fixture in
  `test_stream_extract_annotation_parquets_task.py`.

Run the existing tests for those modules; they should pass unchanged.

- [ ] **Step 1: Write the failing test** (add to the existing test file)

```python
def _ldscore_frame(chrom: int) -> pl.DataFrame:
    return pl.DataFrame(
        {
            POLYFUN_CHR_COL: [chrom, chrom],
            POLYFUN_BP_COL: [1, 2],
            POLYFUN_SNP_COL: ["rsA", "rsB"],
            POLYFUN_A1_COL: ["A", "G"],
            POLYFUN_A2_COL: ["C", "T"],
            "Coding_UCSC_common": [0.1, 0.2],
        }
    )


def test_extracts_only_ldscore_members(tmp_path: Path):
    tarball = tmp_path / "bundle.tar.gz"
    _build_tarball(
        tarball,
        {
            "UKBB_LD/baselineLF2.2.UKB.1.annot.parquet": _member_bytes(_annot_frame(1)),
            "UKBB_LD/baselineLF2.2.UKB.1.l2.ldscore.parquet": _member_bytes(
                _ldscore_frame(1)
            ),
            "UKBB_LD/baselineLF2.2.UKB.2.l2.ldscore.parquet": _member_bytes(
                _ldscore_frame(2)
            ),
            "UKBB_LD/baselineLF2.2.UKB.1.l2.M": b"1\t2\t3\n",
        },
    )
    task = StreamExtractAnnotationParquetsTask(
        meta=_meta(),
        url="http://example.invalid/bundle.tar.gz",
        stream_opener=lambda _url: open(tarball, "rb"),
        required_chromosomes=frozenset({1, 2}),
        member_pattern=LDSCORE_PARQUET_MEMBER_RE,
        dest_stem="baselineLF2.2.UKB.{chrom}.l2.ldscore",
    )
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=_unused_fetch, wf=make_wf())
    names = sorted(p.name for p in result.path.iterdir())
    assert names == [
        "baselineLF2.2.UKB.1.l2.ldscore.parquet",
        "baselineLF2.2.UKB.2.l2.ldscore.parquet",
    ]
```

Import `LDSCORE_PARQUET_MEMBER_RE` from the task module (added in Step 3) and the
`POLYFUN_*_COL` constants from `mecfs_bio.constants.polyfun_constants`.

- [ ] **Step 2: Run test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_stream_extract_annotation_parquets_task.py::test_extracts_only_ldscore_members -v`
Expected: FAIL (`member_pattern`/`dest_stem` are not fields yet;
`LDSCORE_PARQUET_MEMBER_RE` undefined).

- [ ] **Step 3: Generalize the task**

In `stream_extract_annotation_parquets_task.py`: add a module constant and two
`@frozen` fields with defaults preserving current behavior, and use them in
`execute`:

```python
ANNOT_PARQUET_MEMBER_RE = re.compile(r"baselineLF2\.2\.UKB\.(\d+)\.annot\.parquet$")
LDSCORE_PARQUET_MEMBER_RE = re.compile(
    r"baselineLF2\.2\.UKB\.(\d+)\.l2\.ldscore\.parquet$"
)


@frozen
class StreamExtractAnnotationParquetsTask(Task):
    meta: Meta
    url: str
    stream_opener: StreamOpener = field(default=_default_stream_opener, eq=False)
    required_chromosomes: frozenset[int] = frozenset(range(1, 23))
    member_pattern: re.Pattern[str] = ANNOT_PARQUET_MEMBER_RE
    # `{chrom}` is filled with the matched chromosome int; ".parquet" is appended.
    dest_stem: str = "baselineLF2.2.UKB.{chrom}.annot"
```

In `execute`, replace the hardcoded regex/dest with
`self.member_pattern.search(...)` and
`dest = scratch_dir / (self.dest_stem.format(chrom=chrom) + ".parquet")`. Keep the
`.search`, `found[chrom]`, and missing-chromosome logic unchanged. Update the
module docstring to say the task extracts a *selected* member kind (annotation or
LD-score) per the injected pattern — describe current behavior only, no history.

- [ ] **Step 4: Run both stream-extract tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_stream_extract_annotation_parquets_task.py -v`
Expected: PASS (the pre-existing annot test still passes on the defaults).

- [ ] **Step 5: Add the LD-score members asset**

Create `baseline_lf_ldscores.py`:

```python
"""Reference asset: the baseline-LF 2.2.UKB per-chromosome LD-score parquets.

Streamed from the same ~30GB baselineLF_v2.2.UKB.polyfun.tar.gz bundle as the
annotation matrix, keeping the .l2.ldscore.parquet members (each: CHR, SNP, BP,
A1, A2 + 187 annotation-named LD-score columns). The members directory is a
path_remap candidate (~29GB, few files, rarely read).
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.task.annotation_weights.stream_extract_annotation_parquets_task import (
    LDSCORE_PARQUET_MEMBER_RE,
    StreamExtractAnnotationParquetsTask,
)

BASELINE_LF_ANNOTATION_LDSCORE_MEMBERS = StreamExtractAnnotationParquetsTask(
    meta=ReferenceDataDirectoryMeta(
        group="polyfun",
        sub_group="annotations",
        sub_folder=PurePath("raw"),
        id=AssetId("baseline_lf_2.2_ukb_ldscore_parquet_members"),
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.polyfun.tar.gz",
    member_pattern=LDSCORE_PARQUET_MEMBER_RE,
    dest_stem="baselineLF2.2.UKB.{chrom}.l2.ldscore",
)
```

- [ ] **Step 6: Commit**

```bash
git add mecfs_bio/constants/polyfun_constants.py \
        mecfs_bio/build_system/task/annotation_weights/build_baseline_lf_annotation_parquet_task.py \
        mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py \
        mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py \
        mecfs_bio/build_system/task/annotation_weights/stream_extract_annotation_parquets_task.py \
        mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_ldscores.py \
        test_mecfs_bio/unit/build_system/task/annotation_weights/test_stream_extract_annotation_parquets_task.py
git commit -m "feat: extract baseline-LF LD-score members for Approach-2 prior"
```

---

### Task 2: chromosome-blocked weighted ridge submodule

**Files:**
- Create: `mecfs_bio/build_system/task/annotation_weights/chromosome_blocked_ridge.py`
- Test: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_chromosome_blocked_ridge.py`

**Interfaces:**
- Produces (pure numpy; `p` = #features):
  - `ChromRidgeBlock` (`@frozen`): `sw: float`, `swx: np.ndarray (p,)`,
    `swxx: np.ndarray (p,p)`, `swxy: np.ndarray (p,)`, `swy: float`, `swyy: float`.
  - `accumulate_block(x: np.ndarray, y: np.ndarray, w: np.ndarray | None = None) -> ChromRidgeBlock`
  - `combine(blocks: Sequence[ChromRidgeBlock]) -> ChromRidgeBlock`
  - `standardized_system(block) -> StandardizedSystem` (`g_std, b_std, mean, sd, mean_y`)
  - `solve(system, alpha: float) -> np.ndarray` (beta_std)
  - `heldout_r2(held: ChromRidgeBlock, beta_std, train_mean, train_sd, train_mean_y) -> float`
  - `select_alpha_loco(blocks_by_chrom: Mapping[int, ChromRidgeBlock], alphas: Sequence[float]) -> AlphaSelection`
    (`alpha: float`, `mean_r2: float`, `r2_per_chrom: dict[int, float]`)
  - `fit(block, alpha: float) -> RidgeFit` (`beta_raw, beta_std, intercept`)

This is a weighted generalization of the private helpers in
`ridge_annotation_weights_task.py` (`_ChromStats`, `_standardized_system`,
`_solve`, `_heldout_r2`, `_select_alpha_loco`, `_combine`); `w=None` (weights of
1) recovers that unweighted case exactly.

- [ ] **Step 1: Write failing tests**

```python
import numpy as np
import pytest

from mecfs_bio.build_system.task.annotation_weights.chromosome_blocked_ridge import (
    accumulate_block,
    combine,
    fit,
    select_alpha_loco,
)


def _chrom_blocks(seed: int, beta: np.ndarray, n_chrom: int = 6, n_per: int = 200):
    rng = np.random.default_rng(seed)
    blocks = {}
    for chrom in range(1, n_chrom + 1):
        x = rng.normal(size=(n_per, beta.size))
        y = x @ beta + rng.normal(scale=0.01, size=n_per)
        blocks[chrom] = accumulate_block(x, y)
    return blocks


def test_fit_recovers_known_beta_unweighted():
    beta = np.array([2.0, -1.0, 0.5])
    blocks = _chrom_blocks(0, beta)
    result = fit(combine(list(blocks.values())), alpha=1e-6)
    assert np.allclose(result.beta_raw, beta, atol=1e-2)


def test_combine_is_additive():
    beta = np.array([1.0, 0.0])
    blocks = list(_chrom_blocks(1, beta, n_chrom=2).values())
    merged = combine(blocks)
    assert merged.sw == pytest.approx(blocks[0].sw + blocks[1].sw)
    assert np.allclose(merged.swxx, blocks[0].swxx + blocks[1].swxx)


def test_select_alpha_loco_prefers_small_alpha_on_clean_signal():
    beta = np.array([3.0, -2.0, 1.0])
    blocks = _chrom_blocks(2, beta)
    sel = select_alpha_loco(blocks, alphas=(1e-6, 1.0, 1e3, 1e6))
    assert sel.alpha == 1e-6
    assert sel.mean_r2 > 0.999


def test_weights_change_the_fit():
    # A weighted point set whose weighted LS solution differs from unweighted.
    rng = np.random.default_rng(3)
    x = rng.normal(size=(500, 2))
    beta = np.array([1.0, -1.0])
    y = x @ beta + rng.normal(scale=0.1, size=500)
    y[:250] += 5.0  # corrupt half
    w = np.ones(500)
    w[:250] = 1e-6  # downweight the corrupted half
    weighted = fit(accumulate_block(x, y, w), alpha=1e-6)
    unweighted = fit(accumulate_block(x, y), alpha=1e-6)
    assert np.linalg.norm(weighted.beta_raw - beta) < np.linalg.norm(
        unweighted.beta_raw - beta
    )
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_chromosome_blocked_ridge.py -v`
Expected: FAIL (module does not exist).

- [ ] **Step 3: Implement the submodule**

```python
"""Weighted ridge regression accumulated per chromosome block.

Pure numpy: no polars, parquet, or Task imports. A ChromRidgeBlock holds the
weighted cross-product sufficient statistics of one chromosome; combining blocks
lets any leave-one-chromosome-out (or odd/even) fit be formed as a sum, so the
full design matrix is never held in memory. Ridge is fit on centered, per-column
standardized features and unstandardized back; weights of 1 recover ordinary
ridge. See ridge_annotation_weights_task for the unweighted derivations this
generalizes.
"""

from collections.abc import Mapping, Sequence

import numpy as np
from attrs import frozen


@frozen
class ChromRidgeBlock:
    """Weighted cross-product sufficient statistics for one chromosome.

    For the SNPs i in one chromosome, let x_i be the length-p feature row
    (annotation LD-scores), y_i the target (chi-square), and w_i >= 0 the
    regression weight. Every field is a weighted sum over those SNPs:

        sw   = sum_i w_i                     scalar (total weight; = n when w == 1)
        swx  = sum_i w_i x_i                 length p   (weighted feature sum)
        swxx = sum_i w_i x_i x_i^T           p x p      (weighted Gram, X^T W X)
        swxy = sum_i w_i x_i y_i             length p   (weighted feature-target, X^T W y)
        swy  = sum_i w_i y_i                 scalar
        swyy = sum_i w_i y_i^2               scalar

    These are exactly the statistics the weighted-ridge normal equations need,
    and every one is additive across chromosomes, so any leave-one-out or
    odd/even training set is formed by summing the relevant blocks (combine).
    """

    sw: float
    swx: np.ndarray
    swxx: np.ndarray
    swxy: np.ndarray
    swy: float
    swyy: float

    def __attrs_post_init__(self) -> None:
        p = self.swx.shape[0]
        assert self.swx.shape == (p,)
        assert self.swxx.shape == (p, p)
        assert self.swxy.shape == (p,)
        assert self.swx.dtype.kind == "f"


@frozen
class StandardizedSystem:
    """The weighted ridge system after centering and per-column standardization.

    Built from a ChromRidgeBlock. With weighted moments mean_j = swx_j / sw,
    mean_y = swy / sw, and sd_j = sqrt(diag(swxx)_j / sw - mean_j^2) (zeros
    replaced by 1), define the standardized feature z_ij = (x_ij - mean_j)/sd_j.
    Then:

        mean   length p   weighted per-feature mean, mean_j = (sum_i w_i x_ij)/(sum_i w_i)
        sd     length p   weighted per-feature std dev (used to standardize/unstandardize)
        mean_y scalar     weighted mean of the target y
        g_std  p x p      standardized weighted Gram of centered features:
                          g_std[j,k] = (sum_i w_i (x_ij-mean_j)(x_ik-mean_k)) / (sd_j sd_k)
                                     = (swxx - sw * outer(mean, mean))[j,k] / (sd_j sd_k)
        b_std  length p   standardized weighted feature-target cross term:
                          b_std[j] = (sum_i w_i (x_ij-mean_j)(y_i-mean_y)) / sd_j
                                   = (swxy - mean * swy)[j] / sd_j
                          (centering on y is automatic since sum_i w_i (x_ij-mean_j) = 0)

    solve() returns beta_std from (g_std + alpha I) beta_std = b_std; fit() then
    unstandardizes: beta_raw = beta_std / sd, intercept = mean_y - beta_raw . mean.
    """

    g_std: np.ndarray
    b_std: np.ndarray
    mean: np.ndarray
    sd: np.ndarray
    mean_y: float


@frozen
class RidgeFit:
    beta_raw: np.ndarray
    beta_std: np.ndarray
    intercept: float


@frozen
class AlphaSelection:
    alpha: float
    mean_r2: float
    r2_per_chrom: dict[int, float]


def accumulate_block(
    x: np.ndarray, y: np.ndarray, w: np.ndarray | None = None
) -> ChromRidgeBlock:
    x = x.astype(np.float64)
    y = y.astype(np.float64)
    w = np.ones(x.shape[0]) if w is None else w.astype(np.float64)
    wx = x * w[:, None]
    return ChromRidgeBlock(
        sw=float(w.sum()),
        swx=wx.sum(0),
        swxx=x.T @ wx,
        swxy=wx.T @ y,
        swy=float((w * y).sum()),
        swyy=float((w * y * y).sum()),
    )


def combine(blocks: Sequence[ChromRidgeBlock]) -> ChromRidgeBlock:
    return ChromRidgeBlock(
        sw=sum(b.sw for b in blocks),
        swx=sum((b.swx for b in blocks), start=np.zeros_like(blocks[0].swx)),
        swxx=sum((b.swxx for b in blocks), start=np.zeros_like(blocks[0].swxx)),
        swxy=sum((b.swxy for b in blocks), start=np.zeros_like(blocks[0].swxy)),
        swy=sum(b.swy for b in blocks),
        swyy=sum(b.swyy for b in blocks),
    )


def standardized_system(block: ChromRidgeBlock) -> StandardizedSystem:
    n = block.sw
    mean = block.swx / n
    var = np.diag(block.swxx) / n - mean**2
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[sd == 0] = 1.0
    centered_gram = block.swxx - n * np.outer(mean, mean)
    g_std = centered_gram / np.outer(sd, sd)
    b_std = (block.swxy - mean * block.swy) / sd
    return StandardizedSystem(
        g_std=g_std, b_std=b_std, mean=mean, sd=sd, mean_y=block.swy / n
    )


def solve(system: StandardizedSystem, alpha: float) -> np.ndarray:
    p = system.g_std.shape[0]
    return np.linalg.solve(system.g_std + alpha * np.eye(p), system.b_std)


def heldout_r2(
    held: ChromRidgeBlock,
    beta_std: np.ndarray,
    train_mean: np.ndarray,
    train_sd: np.ndarray,
    train_mean_y: float,
) -> float:
    n = held.sw
    c = train_mean_y
    tm, ts = train_mean, train_sd
    ss_res_y = held.swyy - 2.0 * c * held.swy + n * c * c
    z_r = (held.swxy - c * held.swx - tm * held.swy + n * c * tm) / ts
    zz = (
        held.swxx - np.outer(tm, held.swx) - np.outer(held.swx, tm) + n * np.outer(tm, tm)
    ) / np.outer(ts, ts)
    ss_res = ss_res_y - 2.0 * float(beta_std @ z_r) + float(beta_std @ zz @ beta_std)
    mean_y = held.swy / n
    ss_tot = held.swyy - n * mean_y * mean_y
    return 1.0 - ss_res / ss_tot


def select_alpha_loco(
    blocks_by_chrom: Mapping[int, ChromRidgeBlock], alphas: Sequence[float]
) -> AlphaSelection:
    chroms = sorted(blocks_by_chrom)
    best = AlphaSelection(alpha=alphas[0], mean_r2=-np.inf, r2_per_chrom={})
    for alpha in alphas:
        r2s: dict[int, float] = {}
        for held in chroms:
            train = combine([blocks_by_chrom[c] for c in chroms if c != held])
            system = standardized_system(train)
            beta_std = solve(system, alpha)
            r2s[held] = heldout_r2(
                blocks_by_chrom[held],
                beta_std=beta_std,
                train_mean=system.mean,
                train_sd=system.sd,
                train_mean_y=system.mean_y,
            )
        mean_r2 = float(np.mean(list(r2s.values())))
        if mean_r2 > best.mean_r2:
            best = AlphaSelection(alpha=alpha, mean_r2=mean_r2, r2_per_chrom=r2s)
    return best


def fit(block: ChromRidgeBlock, alpha: float) -> RidgeFit:
    system = standardized_system(block)
    beta_std = solve(system, alpha)
    beta_raw = beta_std / system.sd
    intercept = float(system.mean_y - beta_raw @ system.mean)
    return RidgeFit(beta_raw=beta_raw, beta_std=beta_std, intercept=intercept)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_chromosome_blocked_ridge.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/build_system/task/annotation_weights/chromosome_blocked_ridge.py \
        test_mecfs_bio/unit/build_system/task/annotation_weights/test_chromosome_blocked_ridge.py
git commit -m "feat: pure-numpy chromosome-blocked weighted ridge submodule"
```

---

### Task 3: L2-regularized S-LDSC snpvar estimator Task

**Files:**
- Create: `mecfs_bio/build_system/task/annotation_weights/l2_sldsc_snpvar_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_l2_sldsc_snpvar_task.py`

**Interfaces:**
- Consumes: `chromosome_blocked_ridge` (Task 2);
  `BASELINE_LF_ANNOTATION_LDSCORE_MEMBERS` (Task 1);
  `BASELINE_LF_ANNOTATION_MATRIX`; a gwaslab-format REF-oriented build-37 sumstats
  Task carrying `CHR, POS, NEA, EA, Z` (+ optional `N`).
- Produces: `L2RegularizedSldscSnpvarTask` with `.create(asset_id, *,
  sumstats_task, effective_sample_size, annotation_ldscore_members_task=...,
  annotation_matrix_task=BASELINE_LF_ANNOTATION_MATRIX,
  regression_snp_restriction_task=None)`. Output `DirectoryAsset`:
  `snpvar.parquet` (`CHR, POS, NEA, EA, snpvar`) + `diagnostics.json`. Module
  constants: `SNPVAR_PARQUET_FILENAME`, `DIAGNOSTICS_JSON_FILENAME`, `SNPVAR_COL`,
  `MAFBIN_LDSCORE_RE`; private column constants `_CHI2_COL` ("chi2") and
  `_TOTAL_LDSCORE_COL` ("l").

**Algorithm (from spec Component C, with the schema-spike corrections):** total
LD-score `l_i` = sum of the 20 `MAFbin_*` LD-score columns; `M` = reference SNP
count; `chi2 = Z**2` filtered `< 80`; weights
`omega_i = 1 / (het_i * oc_i)`, `het_i = (1 + N * h2bar_parity * l_i / M)**2`,
`oc_i = max(l_i, 1)`; per-parity nested `h2bar`, lambda, and `tau`;
`tau = beta_raw / Nbar`; even chromosomes scored with `tau_odd`, odd with
`tau_even`; `snpvar_i = a_i . tau[opposite_parity]`.

- [ ] **Step 1: Write the failing synthetic test**

Build a tiny synthetic fixture: `n_chrom` chromosomes; two annotations, one
carrying all the heritability. LD-scores `l(i,c)` and annotation values `a(i,c)`
are the same small matrix for the fixture. Simulate `chi2 = 1 + N * (l @ tau) +
noise`, then `Z = sqrt(chi2)`. Assert the recovered `snpvar` ranks variants by
the enriched annotation (Spearman near 1), scores each chromosome from the
opposite parity, and covers 100% of annotation variants.

```python
import numpy as np
import polars as pl
import pytest
from scipy.stats import spearmanr

from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
# ... (Meta/FakeTask/read-spec imports as in test_ridge_annotation_weights_task)
from mecfs_bio.build_system.task.annotation_weights.l2_sldsc_snpvar_task import (
    SNPVAR_COL,
    SNPVAR_PARQUET_FILENAME,
    L2RegularizedSldscSnpvarTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_Z_COL,
)
from mecfs_bio.constants.polyfun_constants import (
    POLYFUN_A1_COL,
    POLYFUN_A2_COL,
    POLYFUN_BP_COL,
    POLYFUN_CHR_COL,
    POLYFUN_SNP_COL,
    POLYFUN_TO_GWASLAB_KEY_RENAME,
)

ANNOTS = ["Coding_UCSC_common", "Repressed_Hoffman_common"]
MAFBINS = [f"MAFbin_lowfreq_{i}" for i in range(1, 11)] + [
    f"MAFbin_frequent_{i}" for i in range(1, 11)
]


def _fixture(tmp_path, n_chrom=6, n_per=300, effective_n=100_000):
    rng = np.random.default_rng(0)
    tau = np.array([2e-7, 0.0])  # only the first annotation carries h2
    ld_rows, ann_rows, ss_rows = [], [], []
    for chrom in range(1, n_chrom + 1):
        a = rng.uniform(0, 1, size=(n_per, len(ANNOTS)))
        # one MAFbin per SNP so the 20 MAFbin ld-scores sum to the total ld-score
        mafbin = np.zeros((n_per, len(MAFBINS)))
        which = rng.integers(0, len(MAFBINS), size=n_per)
        mafbin[np.arange(n_per), which] = rng.uniform(1, 5, size=n_per)
        total_ld = mafbin.sum(1)
        chi2 = 1.0 + effective_n * (a @ tau) + rng.normal(scale=0.05, size=n_per)
        z = np.sqrt(np.maximum(chi2, 1e-6))
        for i in range(n_per):
            key = {POLYFUN_CHR_COL: chrom, POLYFUN_BP_COL: i + 1,
                   POLYFUN_A1_COL: "A", POLYFUN_A2_COL: "G"}
            ld_rows.append({**key, POLYFUN_SNP_COL: f"rs{chrom}_{i}",
                            **{ANNOTS[k]: a[i, k] for k in range(len(ANNOTS))},
                            **{MAFBINS[k]: mafbin[i, k] for k in range(len(MAFBINS))}})
            ann_rows.append({**key, POLYFUN_SNP_COL: f"rs{chrom}_{i}",
                             **{ANNOTS[k]: a[i, k] for k in range(len(ANNOTS))}})
            ss_rows.append({GWASLAB_CHROM_COL: chrom, GWASLAB_POS_COL: i + 1,
                            GWASLAB_NON_EFFECT_ALLELE_COL: "A",
                            GWASLAB_EFFECT_ALLELE_COL: "G", GWASLAB_Z_COL: z[i]})
    # write one ld-score member per chromosome, one annotation matrix, one sumstats
    # (see _make_inputs in test_ridge_annotation_weights_task for the FakeTask +
    # fetch closure pattern). Return paths, tau, ANNOTS.
    ...
```

The assertions:

```python
def test_recovers_enriched_annotation_ranking(tmp_path):
    task, fetch, annots = _fixture(tmp_path)  # builds FakeTasks + fetch closure
    scratch = tmp_path / "scratch"; scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    snpvar = pl.read_parquet(result.path / SNPVAR_PARQUET_FILENAME)
    # the annotation matrix, keyed PolyFun-style; rename to the gwaslab key to join
    ann = pl.read_parquet(annots).rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
    joined = snpvar.join(
        ann, on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL,
                 GWASLAB_NON_EFFECT_ALLELE_COL, GWASLAB_EFFECT_ALLELE_COL]
    )
    rho, _ = spearmanr(joined[SNPVAR_COL], joined["Coding_UCSC_common"])
    assert rho > 0.9
    assert snpvar.height == ann.height  # 100% coverage
```

Add a no-leakage test: perturb the `Z` of one even chromosome, rebuild, and assert
`tau_even` (read from `diagnostics.json`, which stores per-parity tau) is
unchanged while `tau_odd` moves — because even chromosomes only influence the
odd-scoring... (careful with direction: perturbing an even chromosome may change
`tau_even` (fit on even) but must NOT change the snpvar of even chromosomes, which
are scored by `tau_odd`). Concretely assert: perturbing an even chromosome's `Z`
leaves every odd-parity fit output (`tau_even`, the tau used to score even
chromosomes) — i.e. the diagnostics' `tau_odd` — unchanged.

```python
def test_no_leakage_even_perturbation_does_not_move_odd_fit(tmp_path):
    # tau_odd (fit on odd) scores even chromosomes; perturbing an EVEN chromosome
    # must not change tau_odd.
    base = _run_and_read_tau(tmp_path / "a", perturb_even_chrom=None)
    perturbed = _run_and_read_tau(tmp_path / "b", perturb_even_chrom=2)
    assert np.allclose(base["tau_odd"], perturbed["tau_odd"])
    assert not np.allclose(base["tau_even"], perturbed["tau_even"])
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_l2_sldsc_snpvar_task.py -v`
Expected: FAIL (module/class missing).

- [ ] **Step 3: Implement the estimator**

Structure the module with free helpers (per repo convention) and a thin Task.
Key pieces (concrete code for the non-obvious parts):

```python
MAFBIN_LDSCORE_RE = re.compile(r"^MAFbin_(lowfreq|frequent)_\d+$")
SNPVAR_PARQUET_FILENAME = "snpvar.parquet"
DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"
SNPVAR_COL = "snpvar"
# chi2_i = Z_i^2, the regression response.
_CHI2_COL = "chi2"
# l_i = sum over the 20 MAFbin_* LD-score columns = total LD-score of variant i.
_TOTAL_LDSCORE_COL = "l"
_CHI2_CAP = 80.0
_JOIN_KEYS = [
    GWASLAB_CHROM_COL, GWASLAB_POS_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL, GWASLAB_EFFECT_ALLELE_COL,
]


def _ldscore_regression_frame(
    ldscore_path: Path, sumstats: pl.DataFrame
) -> pl.DataFrame:
    """One chromosome's join of LD-scores to sumstats on the four-part key,
    renaming the member's REF-oriented A1/A2 to NEA/EA. Keeps every LD-score
    column (callers select the annotation columns) and adds chi2 and the
    MAFbin-sum total LD-score l."""
    ld = pl.read_parquet(ldscore_path).rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
    mafbins = [c for c in ld.columns if MAFBIN_LDSCORE_RE.match(c)]
    assert len(mafbins) == 20, f"expected 20 MAFbin columns, got {len(mafbins)}"
    ld = ld.with_columns(pl.sum_horizontal(mafbins).alias(_TOTAL_LDSCORE_COL))
    frame = ld.join(sumstats, on=_JOIN_KEYS, how="inner")
    return frame.with_columns(
        (pl.col(GWASLAB_Z_COL) ** 2).alias(_CHI2_COL)
    ).filter(pl.col(_CHI2_COL) < _CHI2_CAP)
```

`execute` orchestration:
1. Load sumstats (`scan_dataframe_asset`) selecting the four keys + `Z` (+ `N` if
   present); drop degenerate `Z` (reuse the repo degenerate-Z guard); collect to a
   polars frame. `Nbar = effective_sample_size` (or `mean(N)` when a per-variant
   `N` column is present).
2. Read `annot_cols` from the annotation-matrix schema (the 187 non-key columns).
   Assert the LD-score members' annotation columns equal `annot_cols` (spec
   decision 6).
3. `M` = total reference SNP count (accumulate row counts while streaming the
   LD-score members in pass A).
4. **Pass A (per-parity h2bar):** for each chromosome member, build the
   regression frame, accumulate univariate sufficient stats of `(_TOTAL_LDSCORE_COL, _CHI2_COL)`
   per chromosome; combine odd / even; `h2bar_parity = max(slope * M / Nbar, 1e-8)`
   where `slope = cov(l, chi2) / var(l)` from the combined univariate sums.
5. **Pass B (per-parity weighted blocks):** for each chromosome member, build the
   regression frame, compute `omega` with `h2bar_parity(chrom)`, and
   `accumulate_block(x=frame[annot_cols], y=frame[_CHI2_COL], w=omega)` into the odd
   or even block dict (keyed by chromosome).
6. `sel_odd = select_alpha_loco(odd_blocks, alphas)`;
   `tau_odd = fit(combine(odd_blocks.values()), sel_odd.alpha).beta_raw / Nbar`;
   symmetric for even.
7. **Pass C (score dense):** stream the annotation matrix per chromosome; for an
   even chromosome use `tau_odd`, for an odd chromosome use `tau_even`;
   `snpvar = annotation_values @ tau`; write `snpvar.parquet` via a streaming
   pyarrow writer with columns `CHR, POS, NEA, EA, snpvar` (rename the annotation
   matrix's PolyFun key via `POLYFUN_TO_GWASLAB_KEY_RENAME`).
8. Write `diagnostics.json`: `alpha` per parity, `mean_heldout_r2` per parity,
   `Nbar`, `M`, `h2bar_odd`, `h2bar_even`, and `tau_odd`/`tau_even` as
   annotation->value maps.

`.create()` derives a `ResultDirectoryMeta` from `sumstats_task.meta` (trait /
project pulled from the dep, per repo convention). Use `alphas` matching
`RidgeAnnotationWeightsTask`'s default grid.

Memory: pass A and B each stream the LD-score members one chromosome at a time
(the ~1.5M x 187 member is the peak, ~1.1GB); the per-chromosome blocks are 22 x
(187x187 f64). Pass C streams the annotation matrix one chromosome at a time.

- [ ] **Step 4: Run tests to verify they pass**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_l2_sldsc_snpvar_task.py -v`
Expected: PASS.

- [ ] **Step 5: Run the wider suite + green**

Run: `pixi r invoke green 2>&1 | tee /tmp/green_task3.log; tail -5 /tmp/green_task3.log`
Expected: no failures (a "no tests ran" from testmon still needs exit 0 AND that
the new tests ran in Step 4).

- [ ] **Step 6: Commit**

```bash
git add mecfs_bio/build_system/task/annotation_weights/l2_sldsc_snpvar_task.py \
        test_mecfs_bio/unit/build_system/task/annotation_weights/test_l2_sldsc_snpvar_task.py
git commit -m "feat: L2-regularized S-LDSC snpvar estimator (Approach-2 prior)"
```

---

### Task 4: wire the Approach-2 prior into the fine-mapping generator

**Files:**
- Modify: `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py`
- Create: a trait wiring instance (e.g. DecodeME) under
  `mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/` — mirror the existing
  Approach-1 fine-mapping wiring for that trait.

**Interfaces:**
- Produces: a `@frozen PolyfunPriorSource(prior_task: Task, prior_col: str,
  prior_pipe: DataProcessingPipe)` with a module default built from
  `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS`, `POLYFUN_PRIOR_COL`,
  `create_prior_col_pipe(q_factor)`. `generate_assets_polyfun_explain_fine_map`
  and `generate_polyfun_explain_group` accept `prior_source: PolyfunPriorSource`
  (default preserves current behavior byte-for-byte).

Per repo convention there are **no structural tests for asset-generator wiring**
(import-time construction suffices); correctness is covered by Task 3's tests and
the coverage assertion.

- [ ] **Step 1: Add the `PolyfunPriorSource` value object and thread it through**

In the generator, replace the hardcoded `PriorInfo(...)` construction in
`generate_polyfun_explain_group` with one built from `prior_source`:

```python
@frozen
class PolyfunPriorSource:
    prior_task: Task
    prior_col: str
    prior_pipe: DataProcessingPipe


def default_polyfun_prior_source(q_factor: int) -> PolyfunPriorSource:
    return PolyfunPriorSource(
        prior_task=COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
        prior_col=POLYFUN_PRIOR_COL,
        prior_pipe=create_prior_col_pipe(q_factor),
    )
```

`SharedFineMapInputs` gains a `prior_source: PolyfunPriorSource` field;
`generate_polyfun_explain_group` builds
`PriorInfo(prior_task=shared.prior_source.prior_task,
prior_col=shared.prior_source.prior_col,
prior_pipe=shared.prior_source.prior_pipe)`.
`generate_assets_polyfun_explain_fine_map` gains
`prior_source: PolyfunPriorSource | None = None`, defaulting via
`default_polyfun_prior_source(q_factor)` when `None`.

- [ ] **Step 2: Verify the default path is unchanged**

Run: `pixi r pytest test_mecfs_bio -k polyfun_explain -q`
Expected: PASS (unchanged behavior). Also `pixi r python -c "import
mecfs_bio.asset_generator.polyfun_explain_fine_mapping_asset_generator"` imports
clean.

- [ ] **Step 3: Add the Approach-2 wiring instance for one trait**

Construct an `L2RegularizedSldscSnpvarTask` for the trait's build-37 sumstats,
and a `PolyfunPriorSource(prior_task=that_task, prior_col=SNPVAR_COL,
prior_pipe=create_prior_col_pipe(q_factor))`, and call
`generate_assets_polyfun_explain_fine_map(..., prior_source=approach2_source)` for
the same locus already fine-mapped under Approach 1, so their credible sets are
comparable via the existing UpSet-over-runs tooling.

- [ ] **Step 4: green + commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_task4.log; tail -5 /tmp/green_task4.log
git add mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py \
        mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/
git commit -m "feat: injectable prior source; wire Approach-2 prior into fine-mapping"
```

---

### Task 5 (separate): migrate `RidgeAnnotationWeightsTask` onto the submodule

Independent of Tasks 1-4; do it only after they are green, so the new feature
never carries the refactor's risk. Gated by a numeric-equivalence test.

**Files:**
- Modify: `mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py` (add equivalence test)

- [ ] **Step 1: Add a numeric-equivalence test**

Capture the current output on a rigged system BEFORE refactoring: the existing
`test_recovers_known_linear_weights` already pins `gamma_raw ~= truth` and
`mean_heldout_r2 > 0.999`. Add an explicit regression test that records the fitted
`gamma_raw` vector for a fixed seed and asserts it is reproduced to `atol=1e-10`
after the refactor (compute the expected vector by running the current code once
and pasting the values, or assert equality against a re-run through the submodule
on the same rigged blocks).

- [ ] **Step 2: Refactor to delegate to the submodule**

Replace `_ChromStats`, `_standardized_system`, `_solve`, `_heldout_r2`,
`_select_alpha_loco`, `_combine` with calls into `chromosome_blocked_ridge`
(`accumulate_block(x, y)` with no weights, `select_alpha_loco`, `fit`). Keep
`_accumulate_per_chromosome` producing `ChromRidgeBlock`s per chromosome; keep the
task's output schema (`weights.parquet`, `diagnostics.json`) identical.

- [ ] **Step 3: Run the equivalence + existing tests**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py -v`
Expected: PASS (bit-identical coefficients).

- [ ] **Step 4: green + commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/green_task5.log; tail -5 /tmp/green_task5.log
git add mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py \
        test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py
git commit -m "refactor: RidgeAnnotationWeightsTask onto shared chromosome-blocked ridge"
```

---

## Self-Review

**Spec coverage:** Component A → Task 1; Component B → Task 2; Component C → Task
3; Component D → Task 4; the gated `RidgeAnnotationWeightsTask` migration → Task 5.
Decisions 1-8: pragmatic weighted ridge (Task 2/3); dense streamed regression
(Task 3 passes); shared submodule (Task 2, migration Task 5); MAFbin-sum
`l_i`/`M` (Task 3 `_ldscore_regression_frame`); floor at the seam (Task 4
`create_prior_col_pipe`); column-name match assert (Task 3 Step 3.2); exact
four-part join key (global constraint, `_JOIN_KEYS`); nested per-parity
h2bar/lambda/tau + no-leakage test (Task 3). Out-of-scope items (Approach 3,
custom annotations, tau-based explanation) are not planned, as intended.

**Placeholder scan:** the estimator `execute` orchestration (Task 3 Step 3) is
given as an ordered spec of concrete operations plus concrete code for the
non-obvious join/weight helper; the fixture builder is sketched with the
`...`-marked file-writing tail that mirrors `_make_inputs` in the existing ridge
test — the executor completes it against that pattern. No "TBD"/"handle edge
cases"/"similar to Task N".

**Type consistency:** `ChromRidgeBlock`/`RidgeFit`/`AlphaSelection` fields and the
`accumulate_block`/`combine`/`fit`/`select_alpha_loco` signatures used in Task 3
match Task 2's definitions; `SNPVAR_COL`, `SNPVAR_PARQUET_FILENAME` referenced in
Task 4 match Task 3's constants; the four-part `_JOIN_KEYS` uses the gwaslab
constants named in Global Constraints.

## Execution note (worktree)

A background run may be using the main git dir. If so, create an isolated
workspace with superpowers:using-git-worktrees before executing; otherwise this
branch (`polyfun-approach-2-l2-sldsc-prior`) in the main dir is fine (confirm with
the user).
