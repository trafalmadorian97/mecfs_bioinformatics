# Polyfun-explain annotation callouts — implementation plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add per-credible-set on-figure callouts to the polyfun-explain plot that name, for the variants the prior most clearly boosted, the annotation families driving the boost — e.g. `173845678:A:T (cons ++, cod +)`.

**Architecture:** All selection and label text is computed in `PolyfunExplainContrastTask` (deterministic, table-out) and written to a new `callouts.parquet`. `PolyfunExplainPlotTask` reads that table and only places the labels with `textalloc`, deriving no statistics of its own.

**Tech Stack:** Python, polars, numpy, matplotlib (OO API), textalloc, pixi (`pixi r ...`).

**Spec:** `experiments/claude/design_specs/2026-08-30-polyfun-explain-annotation-callouts-design.md`

## Global Constraints

- Run everything via pixi: `pixi r python ...`, `pixi r invoke green`.
- Light testing only: one small unit test on the pure label helpers; rely on the demo rebuild for end-to-end. Do not add heavy/rigged fixtures.
- Repo conventions: column-name constants (no repeated string literals), docstrings without backticks/RST, prefer polars/Path, prefer free helper functions taking only what they need, named kwargs for same-typed params.
- `invoke green` (ruff lintfix + format, ty typecheck, pytest-testmon) must pass at exit 0 before each commit.
- Commit trailers on every commit:
  ```
  Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01AbyeBRGx8SnUcj6YtmpB3m
  ```

---

## File Structure

- **Modify** `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py`
  - New public constants: `CALLOUTS_FILENAME`, `CALLOUT_CS_COL`, `CALLOUT_PIP_PF_COL`, `CALLOUT_PIP_U_COL`, `CALLOUT_LABEL_COL`, and a `_CALLOUT_SCHEMA`.
  - New free helpers: `_family_background_sd`, `_callout_families`, `_format_callout_label`, `_build_callouts`.
  - Wire callout computation into `execute()` and write `callouts.parquet`.
- **Create** `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_callouts.py`
  - Unit tests for the two pure helpers `_callout_families` and `_format_callout_label`.
- **Modify** `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py`
  - Read `callouts.parquet` in `execute()`, pass to `_render`, place labels on the polyfun PIP panel with `textalloc`, add panel headroom.
- **Modify** `experiments/claude/polyfun_explain_plot/rebuild_demo_plot.py` — no change expected; used to rebuild + view.

---

## Task 1: Compute callouts in the contrast task

**Files:**
- Modify: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_callouts.py`

**Interfaces:**
- Consumes (all already present in the module): `_KEY`, `PIP_COLUMN`, `FAMILY_COL`, `FAMILY_CONTRAST_COL`, `FAMILY_SCALED_COL`, `_family_scaled`, `_load_cs_numbers` (yields `_KEY` + `cs_number`), `GWASLAB_POS_COL`, `GWASLAB_EFFECT_ALLELE_COL`, `GWASLAB_NON_EFFECT_ALLELE_COL`, `GWASLAB_CHROM_COL`; from `mecfs_bio.constants.polyfun_annotation_families`: `FAMILY_SHORT_LABELS`, `AnnotationFamily`.
- Produces (later tasks rely on these):
  - `CALLOUTS_FILENAME = "callouts.parquet"`.
  - `callouts.parquet` columns: `_KEY` (CHR:Int64, POS:Int64, EA:str, NEA:str), `CALLOUT_CS_COL="cs"` (Int32), `CALLOUT_PIP_PF_COL="pip_pf"` (Float64), `CALLOUT_PIP_U_COL="pip_u"` (Float64), `CALLOUT_LABEL_COL="label"` (str).
  - `_format_callout_label(focal_key: dict, families: list[tuple[str, str]]) -> str`
  - `_callout_families(per_family: pl.DataFrame, focal_key: dict, family_sd: dict[str, float], max_families: int) -> list[tuple[str, str]]` returning `(family_name, marker)` pairs, marker in `{"+", "++"}`, ordered by z descending.

- [ ] **Step 1: Add constants** near the other filename/column constants (after `SELECTION_IMPORTANT_FAMILIES_KEY`):

```python
CALLOUTS_FILENAME = "callouts.parquet"
CALLOUT_CS_COL = "cs"
CALLOUT_PIP_PF_COL = "pip_pf"
CALLOUT_PIP_U_COL = "pip_u"
CALLOUT_LABEL_COL = "label"
# Fixed schema so an empty callout set still round-trips through parquet and the
# plot task can read a well-typed (possibly zero-row) frame.
_CALLOUT_SCHEMA: dict[str, pl.DataType] = {
    GWASLAB_CHROM_COL: pl.Int64,
    GWASLAB_POS_COL: pl.Int64,
    GWASLAB_EFFECT_ALLELE_COL: pl.String,
    GWASLAB_NON_EFFECT_ALLELE_COL: pl.String,
    CALLOUT_CS_COL: pl.Int32,
    CALLOUT_PIP_PF_COL: pl.Float64,
    CALLOUT_PIP_U_COL: pl.Float64,
    CALLOUT_LABEL_COL: pl.String,
}
# Selection thresholds (see design doc). Change-based only; no absolute PIP floor.
_DOMINANCE_MARGIN = 0.05
_PRIOR_EFFECT_MARGIN = 0.10
_MAX_CALLOUT_FAMILIES = 3
```

- [ ] **Step 2: Add the pure label + family-selection helpers** (place with the other free helpers, e.g. after `_select_families`):

```python
def _format_callout_label(
    focal_key: dict, families: list[tuple[str, str]]
) -> str:
    """Render one callout string: pos:nea:ea, then the key families with their
    strength markers. No families -> just the variant id (no parentheses)."""
    head = (
        f"{focal_key[GWASLAB_POS_COL]}:"
        f"{focal_key[GWASLAB_NON_EFFECT_ALLELE_COL]}:"
        f"{focal_key[GWASLAB_EFFECT_ALLELE_COL]}"
    )
    if not families:
        return head
    inner = ", ".join(
        f"{FAMILY_SHORT_LABELS[fam]} {marker}" for fam, marker in families
    )
    return f"{head} ({inner})"


def _callout_families(
    per_family: pl.DataFrame,
    focal_key: dict,
    family_sd: dict[str, float],
    max_families: int,
) -> list[tuple[str, str]]:
    """The key families for one flagged variant: those whose per-family contrast
    is positive (elevated at this variant) and exceeds one background SD, bucketed
    1-2 SD -> '+', >2 SD -> '++'. Top max_families by z, z descending. Families
    with a degenerate (<= 0) background SD are skipped."""
    focal = per_family
    for k, v in focal_key.items():
        focal = focal.filter(pl.col(k) == v)
    scored: list[tuple[float, str, str]] = []
    for row in focal.iter_rows(named=True):
        fam = row[FAMILY_COL]
        diff = row[FAMILY_CONTRAST_COL]
        sd = family_sd.get(fam, 0.0)
        if diff <= 0.0 or sd <= 0.0:
            continue
        z = diff / sd
        if z <= 1.0:
            continue
        scored.append((z, fam, "++" if z > 2.0 else "+"))
    scored.sort(key=lambda t: t[0], reverse=True)
    return [(fam, marker) for _, fam, marker in scored[:max_families]]
```

- [ ] **Step 3: Add the background-SD helper** (uniform-PIP-weighted SD of `family_scaled` per family over the uniform background — the same population and weights used for `abar`, so the mean it implies equals `sum_c gamma_c*abar_c`):

```python
def _family_background_sd(
    uni_annot: pl.DataFrame,
    annot_cols: list[str],
    gamma: dict[str, float],
    family: dict[str, str],
) -> dict[str, float]:
    """Per family, the uniform-PIP-weighted standard deviation of family_scaled
    over the uniform-run variants. Falls back to equal weights when the uniform
    run carried no signal (total PIP <= 0), matching the abar fallback."""
    fs = _family_scaled(uni_annot, annot_cols, gamma, family)
    weight = uni_annot.select(*_KEY, pl.col(PIP_COLUMN).alias("w"))
    if uni_annot[PIP_COLUMN].sum() <= 0.0:
        weight = weight.with_columns(pl.lit(1.0).alias("w"))
    stats = (
        fs.join(weight, on=_KEY, how="inner")
        .group_by(FAMILY_COL)
        .agg(
            (pl.col("w") * pl.col(FAMILY_SCALED_COL)).sum().alias("wx"),
            (pl.col("w") * pl.col(FAMILY_SCALED_COL) ** 2).sum().alias("wx2"),
            pl.col("w").sum().alias("wsum"),
        )
    )
    out: dict[str, float] = {}
    for row in stats.iter_rows(named=True):
        wsum = row["wsum"]
        if wsum <= 0.0:
            continue
        mean = row["wx"] / wsum
        var = max(row["wx2"] / wsum - mean * mean, 0.0)
        out[row[FAMILY_COL]] = float(np.sqrt(var))
    return out
```

- [ ] **Step 4: Add the callout builder** (iterates polyfun credible sets, applies the two gates, builds a row per flagged variant):

```python
def _build_callouts(
    pf_variants: pl.DataFrame,
    uni_variants: pl.DataFrame,
    cs_pf: pl.DataFrame,
    per_family: pl.DataFrame,
    family_sd: dict[str, float],
) -> pl.DataFrame:
    """One callout row per polyfun credible set whose top-PIP variant clears both
    gates: PIP >= _DOMINANCE_MARGIN above the next-highest PIP in the same CS, and
    PIP >= _PRIOR_EFFECT_MARGIN above the same variant's uniform PIP (0 if absent
    from the uniform run)."""
    cs = cs_pf.join(
        pf_variants.select(*_KEY, PIP_COLUMN), on=_KEY, how="inner"
    )
    uni_pip = {
        tuple(row[k] for k in _KEY): row[PIP_COLUMN]
        for row in uni_variants.select(*_KEY, PIP_COLUMN).iter_rows(named=True)
    }
    rows: list[dict] = []
    for (cs_number,), grp in cs.group_by("cs_number", maintain_order=True):
        grp = grp.sort(PIP_COLUMN, descending=True)
        top = grp.row(0, named=True)
        top_pip = top[PIP_COLUMN]
        next_pip = grp[PIP_COLUMN][1] if grp.height > 1 else 0.0
        if top_pip - next_pip < _DOMINANCE_MARGIN:
            continue
        focal_key = {k: top[k] for k in _KEY}
        u = uni_pip.get(tuple(top[k] for k in _KEY), 0.0)
        if top_pip - u < _PRIOR_EFFECT_MARGIN:
            continue
        families = _callout_families(
            per_family, focal_key, family_sd, _MAX_CALLOUT_FAMILIES
        )
        rows.append(
            {
                **{k: focal_key[k] for k in _KEY},
                CALLOUT_CS_COL: int(cs_number),
                CALLOUT_PIP_PF_COL: float(top_pip),
                CALLOUT_PIP_U_COL: float(u),
                CALLOUT_LABEL_COL: _format_callout_label(focal_key, families),
            }
        )
    return pl.DataFrame(rows, schema=_CALLOUT_SCHEMA)
```

- [ ] **Step 5: Wire into `execute()`.** After `family_scaled = _family_scaled(...)` and the existing focal/display block, before `return DirectoryAsset(scratch_dir)`, add:

```python
        family_sd = _family_background_sd(uni_annot, annot_cols, gamma, family)
        callouts = _build_callouts(
            pf_variants=pf_variants,
            uni_variants=uni_variants,
            cs_pf=cs_pf,
            per_family=per_family,
            family_sd=family_sd,
        )
        callouts.write_parquet(scratch_dir / CALLOUTS_FILENAME)
```

Note: `uni_annot`, `annot_cols`, `gamma`, `family`, `pf_variants`, `uni_variants`, `cs_pf`, `per_family` are all already local in `execute()`. `pf_variants` already carries `PIP_COLUMN` (from `_load_run_variants`) and `DISP_LIFT`.

- [ ] **Step 6: Write the unit test** for the two pure helpers:

```python
import polars as pl

from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    FAMILY_CONTRAST_COL,
    _callout_families,
    _format_callout_label,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    FAMILY_COL,
)

_KEY_COLS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


def _key(pos: int) -> dict:
    return {
        GWASLAB_CHROM_COL: 1,
        GWASLAB_POS_COL: pos,
        GWASLAB_EFFECT_ALLELE_COL: "T",
        GWASLAB_NON_EFFECT_ALLELE_COL: "A",
    }


def test_callout_families_buckets_and_orders_by_z():
    focal = _key(123)
    per_family = pl.DataFrame(
        [
            {**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 3.0},   # z=3 -> ++
            {**focal, FAMILY_COL: "coding", FAMILY_CONTRAST_COL: 1.5},      # z=1.5 -> +
            {**focal, FAMILY_COL: "repressed", FAMILY_CONTRAST_COL: 0.5},   # z=0.5 -> drop
            {**focal, FAMILY_COL: "histone_marks", FAMILY_CONTRAST_COL: -4.0},  # negative -> drop
        ]
    )
    family_sd = {
        "conserved": 1.0,
        "coding": 1.0,
        "repressed": 1.0,
        "histone_marks": 1.0,
    }
    result = _callout_families(per_family, focal, family_sd, max_families=3)
    assert result == [("conserved", "++"), ("coding", "+")]


def test_callout_families_skips_degenerate_sd():
    focal = _key(123)
    per_family = pl.DataFrame(
        [{**focal, FAMILY_COL: "conserved", FAMILY_CONTRAST_COL: 3.0}]
    )
    assert _callout_families(per_family, focal, {"conserved": 0.0}, 3) == []


def test_format_label_with_and_without_families():
    focal = _key(174128548)
    assert (
        _format_callout_label(focal, [("conserved", "++"), ("coding", "+")])
        == "174128548:A:T (cons ++, cod +)"
    )
    assert _format_callout_label(focal, []) == "174128548:A:T"
```

- [ ] **Step 7: Run the test.**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_callouts.py -q`
Expected: 3 passed.

- [ ] **Step 8: Run green, then commit.**

Run: `pixi r invoke green 2>&1 | tee /tmp/green.log | tail -5` (confirm exit 0, ty passes).

```bash
git add mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_callouts.py
git commit -m "feat(polyfun-explain): compute per-CS annotation callouts in contrast task

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01AbyeBRGx8SnUcj6YtmpB3m"
```

---

## Task 2: Render callouts on the polyfun PIP panel

**Files:**
- Modify: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py`

**Interfaces:**
- Consumes from Task 1: `CALLOUTS_FILENAME`, `CALLOUT_PIP_PF_COL`, `CALLOUT_LABEL_COL`, and the callout frame's `GWASLAB_POS_COL` column.
- Produces: the figure now carries labels; no new external interface.

- [ ] **Step 1: Add imports.** Add `import textalloc as ta` with the other third-party imports, and extend the contrast-task import to pull the callout constants:

```python
import textalloc as ta
```

```python
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_contrast_task import (
    CALLOUT_LABEL_COL,
    CALLOUT_PIP_PF_COL,
    CALLOUTS_FILENAME,
    _load_run_variants,
)
```

- [ ] **Step 2: Load the callouts in `execute()`** (the contrast dir is a dep already). After `ld = np.load(...)` add:

```python
        contrast_dir = _dir(fetch, self.contrast_task)
        callouts = pl.read_parquet(contrast_dir / CALLOUTS_FILENAME)
```

and pass `callouts=callouts` into the `_render(...)` call.

- [ ] **Step 3: Thread `callouts` through `_render`.** Add `callouts: pl.DataFrame` to the `_render` signature (e.g. after `pf_full`).

- [ ] **Step 4: Add headroom + label placement.** Replace the shared-PIP-scale block

```python
    pip_top = _shared_pip_top(uni_cs, pf_cs)
    axes[1].set_ylim(0.0, pip_top)
    axes[2].set_ylim(0.0, pip_top)
```

with a version that reserves the top ~40% of both PIP panels for labels (both raised equally so the two panels keep a shared data scale and stay comparable), then places the callout labels on the polyfun panel:

```python
    pip_top = _shared_pip_top(uni_cs, pf_cs)
    # Reserve headroom above the stems for the callout labels; raise BOTH PIP
    # panels equally so their data scale stays shared (directly comparable).
    label_top = pip_top / 0.6
    axes[1].set_ylim(0.0, label_top)
    axes[2].set_ylim(0.0, label_top)
    _place_callouts(axes[2], callouts, bp_min, bp_max)
```

- [ ] **Step 5: Add the `_place_callouts` helper** (bottom of the file):

```python
def _place_callouts(
    ax_pf, callouts: pl.DataFrame, bp_min: int, bp_max: int
) -> None:
    """Annotate the polyfun PIP panel: one text label per callout row, anchored at
    (POS, pip_pf), with textalloc arranging them to avoid mutual overlap and the
    stems, joined to their anchors by thin leader lines. Empty frame -> no-op."""
    if callouts.height == 0:
        return
    xs = callouts[GWASLAB_POS_COL].to_numpy().astype(float).tolist()
    ys = callouts[CALLOUT_PIP_PF_COL].to_numpy().astype(float).tolist()
    texts = callouts[CALLOUT_LABEL_COL].to_list()
    # Seed so textalloc's candidate search is reproducible across builds (keeps
    # the committed SVG stable); textalloc draws from numpy's global RNG.
    np.random.seed(0)
    ta.allocate(
        ax_pf,
        xs,
        ys,
        texts,
        x_scatter=xs,
        y_scatter=ys,
        textsize=7,
        linecolor="black",
        linewidth=0.6,
        avoid_label_lines_overlap=True,
    )
```

- [ ] **Step 6: Update the module docstring** panel list to mention the polyfun PIP panel carries per-CS callouts, and drop the stale "forthcoming per-variant annotation label" note.

- [ ] **Step 7: Run the existing plot smoke test.**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_plot_task.py -q`
Expected: 1 passed. (The synthetic contrast inputs now also write `callouts.parquet`, which the plot reads; an empty frame renders today's figure.)

- [ ] **Step 8: Run green, then commit.**

Run: `pixi r invoke green 2>&1 | tee /tmp/green.log | tail -5` (exit 0, ty passes).

```bash
git add mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py
git commit -m "feat(polyfun-explain): render per-CS annotation callouts on the polyfun PIP panel

Co-Authored-By: Claude Opus 4.8 <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01AbyeBRGx8SnUcj6YtmpB3m"
```

---

## Task 3: Rebuild the demo and inspect

**Files:** none (uses `experiments/claude/polyfun_explain_plot/rebuild_demo_plot.py`).

- [ ] **Step 1: Force-rebuild the l1 demonstrator plot** (contrast + plot are code-changed, so force transitive rebuild):

Run: `pixi r python experiments/claude/polyfun_explain_plot/rebuild_demo_plot.py 2>&1 | tail -5`
Expected: ends with `done: DirectoryAsset(...)` and the plot asset path.

- [ ] **Step 2: View the PNG** at
`assets/base_asset_store/gwas/ME_CFS/DecodeME/analysis/decode_me_polyfun_explainchr1_173500000_174500000_palindromes_keep_l1_explain_plot/explain_plot.png`
and confirm: the polyfun PIP panel shows a callout on its top variant with a `pos:nea:ea (families)` label, leader line to the stem, no overlap with stems or axis. (L=1 has one CS, so expect one callout if both gates pass; on this locus the top polyfun PIP ~0.32 vs uniform ~0.03 and next-in-CS ~0.15, so both gates pass.)

- [ ] **Step 3: Sanity-check an L>1 config** to see multiple callouts / collision handling. Temporarily point the demo at the `l10` group (`POLYFUN_EXPLAIN_CHR1_174.groups[2].plot`), rebuild, view, then revert the demo script. (No commit; inspection only.)

- [ ] **Step 4: Report to the user** with the viewed image and a one-line summary of what the callouts show, and ask whether to keep this plot approach (per the stated goal). Do not proceed to multi-locus rollout or further polish unprompted.

---

## Self-Review

**Spec coverage:**
- Selection gates (dominance 5pp, prior-effect 10pp vs full uniform PIP, no floor) → Task 1 Step 4 (`_build_callouts`, `_DOMINANCE_MARGIN`, `_PRIOR_EFFECT_MARGIN`, `uni_pip` default 0).
- One-per-CS (self-limiting) → Task 1 Step 4 iterates per `cs_number`, takes the single top variant.
- Key families per callout, positive-only, z-bucketed +/++ , top-3 → Task 1 Step 2 (`_callout_families`).
- Background SD (uniform-PIP-weighted, same population as abar), degenerate-SD skip → Task 1 Step 3 (`_family_background_sd`), guarded in `_callout_families`.
- Empty-families label with no parentheses → Task 1 Step 2 (`_format_callout_label`), tested Step 6.
- `callouts.parquet` schema → Task 1 Step 1 (`_CALLOUT_SCHEMA`), written Step 5.
- Compute in contrast task, plot reads it → Task 1 (compute), Task 2 (read + render).
- textalloc on polyfun PIP panel, headroom, determinism → Task 2 Steps 4-5.
- Light testing → one unit-test file (Task 1 Step 6); plot keeps its smoke test.

**Placeholder scan:** none — every step has concrete code or an exact command.

**Type consistency:** `_callout_families` returns `list[tuple[str, str]]` consumed by `_format_callout_label` (Task 1) — matches. `CALLOUTS_FILENAME`, `CALLOUT_PIP_PF_COL`, `CALLOUT_LABEL_COL` defined in Task 1 Step 1, imported in Task 2 Step 1 — names match. `_build_callouts` reads `pf_variants[PIP_COLUMN]` / `cs_pf["cs_number"]` — both confirmed present (`_load_run_variants`, `_load_cs_numbers`).

**Open risk to watch during execution:** `ta.allocate`'s exact keyword set / RNG behavior may differ from the magma call site — if a kwarg is unsupported, drop it; if placement still churns the SVG, that only affects the plot asset hash (nothing consumes it downstream), so it does not block.
