# Retire the Unordered Allele Key Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the SNV-only `unordered_allele_key` with exact `(chrom, pos, ea, nea)` joins on the polyfun/SUSIE fine-mapping path, make the remaining unordered-key uses SNV-safe by assertion, and enable the harmonization build-match guard in SUSIE.

**Architecture:** The unordered allele key silently mis-keys indels (mirrored `T/TCA` vs `TCA/T` collide). Everything on the fine-mapping path is now verified reference-oriented (`NEA == REF`), so exact 4-tuple joins are both safe and strictly more correct. This plan: (1) adds an `assert_all_snv` guard and applies it where the key stays (PPP index) or the task is SNV-only (the allele harmonizer); (2) converts the baseline-LF annotation join path (producer + ridge weights + explain contrast) and the SUSIE polyfun-prior join to exact keys; (3) removes `HarmonizeGWASWithReferenceViaAlleles` from the two fine-mapping generators, moving its palindrome/pipe responsibilities into SUSIE's `gwas_data_pipe`; and (4) enables the construction-time build-match assertion in `SusieRFinemapTask.create`.

**Tech Stack:** Python 3.13, attrs (`@frozen`), polars + narwhals pipes, R/rpy2 (SUSIE), the repo's `Task`/`Asset`/`Meta` build system, pixi, pytest.

**Spec:** `experiments/claude/design_specs/susie_polyfun_unordered_allele_key_cleanup.md`

**Depends on:** the harmonization-provenance plan (`2026-09-18-harmonization-provenance-plumbing-plan.md`) must be merged first — Task 7 here consumes `harmonization_info` that plan threads through the metadata, and Task 6 assumes the LD labels are renamed via `RenameColsTask` (Plan 1 Task 8).

## Global Constraints

- All commands via pixi: `pixi r <cmd>`, `pixi r python <script>`, `pixi r invoke green`.
- After any significant change run `pixi r invoke green`. testmon skips unaffected tests; name a test in `pixi r pytest ... -v` to see it run.
- Column-name constants from `mecfs_bio.constants.gwaslab_constants` (`GWASLAB_CHROM_COL="CHR"`, `GWASLAB_POS_COL="POS"`, `GWASLAB_EFFECT_ALLELE_COL="EA"`, `GWASLAB_NON_EFFECT_ALLELE_COL="NEA"`). Never repeat the string literals.
- Verified orientation fact this plan relies on: annotation `A1` == 1kg-prior `A1` == LD `allele1` == gwas `NEA` == hg19 FASTA REF (see `experiments/claude/orientation_check/verify_fasta_orientation.py`). So the annotation/prior column `A1` maps to `NEA`, `A2` maps to `EA`.
- Prefer helper free functions; Path over str; assert-narrow over cast; no monkeypatching. Test stub is `FakeTask(meta=...)` from `mecfs_bio.build_system.task.fake_task`.
- **All imports at the top of the file** — never inside a function or test body.
- Docstrings: no backticks around inline code, no RST.
- **Never edit text under `docs/`** — it is the human-written published site. If a change alters a documented result, surface the proposed doc update to the user; do not edit `docs/` yourself.
- **Do not run the from-scratch system test locally.** It is slow; run it in GitHub Actions on the PR branch as the final check.

### Testing philosophy
Test only real logic with a real failure mode. New tests are added only for: `assert_all_snv` (the guard), `DropPalindromesPipe` (the filter), and the SUSIE build-match assertion (the guard) — plus updating the existing contrast unit test where its fixtures encode the old allele-key behavior. No field-you-just-set, thin-wrapper, or asset-generator tests.

The result-affecting changes are guarded as follows:
- Annotation exact join (Task 3): the unit test `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py` (run locally).
- Prior exact join (Task 4) and generator harmonizer removal (Task 6): the from-scratch system test `test_mecfs_bio/system/test_decode_me_polyfun_susie.py`, **run in CI on the PR, not locally**. Note it exercises the **keep**-palindrome path only — its documented lead variant is the palindrome chr1:173,855,298 {A,T}. The `drop`-palindrome path is not separately guarded, which is acceptable: `keep` is the primary setting (we already know the input is correctly reference-harmonized, so palindromes are resolvable and normally kept).

### Ordering (dependencies)
Tasks 1→8 in order. Critical: Task 6 (remove harmonizer from generators) MUST precede Task 8 (SNV-assert the harmonizer), because the fine-mapping gwas contains indels and would trip the assert while the generators still feed it. Task 5 (pipes) precedes Task 6. Task 7 (SUSIE build-match) follows Task 6 (only then does `harmonization_info` reach SUSIE).

---

### Task 1: `assert_all_snv` guard + document `unordered_allele_key`

**Files:**
- Modify: `mecfs_bio/build_system/task/ppp_database/allele_key.py`
- Test: `test_mecfs_bio/unit/build_system/task/ppp_database/test_allele_key_snv_guard.py`

**Interfaces:**
- Produces: `assert_all_snv(df: pl.DataFrame, *allele_cols: str) -> None` — raises `AssertionError` if any value in any named column is not a single base.

- [ ] **Step 1: Write the failing test**

```python
# test_mecfs_bio/unit/build_system/task/ppp_database/test_allele_key_snv_guard.py
import polars as pl
import pytest

from mecfs_bio.build_system.task.ppp_database.allele_key import assert_all_snv


def test_assert_all_snv_passes_for_snvs():
    assert_all_snv(pl.DataFrame({"EA": ["A", "C"], "NEA": ["G", "T"]}), "EA", "NEA")


def test_assert_all_snv_rejects_indels():
    df = pl.DataFrame({"EA": ["A", "TCA"], "NEA": ["G", "T"]})
    with pytest.raises(AssertionError):
        assert_all_snv(df, "EA", "NEA")
```

- [ ] **Step 2: Run to confirm failure**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/ppp_database/test_allele_key_snv_guard.py -v`
Expected: FAIL (ImportError).

- [ ] **Step 3: Implement the guard and document the key**

In `allele_key.py`, extend the module docstring / `unordered_allele_key` docstring to state it is valid only for SNVs (mirrored indels such as `T/TCA` and `TCA/T` collide on the sorted key), and add:

```python
def assert_all_snv(df: pl.DataFrame, *allele_cols: str) -> None:
    """Fail fast if any allele in the named columns is not a single base.

    unordered_allele_key and the PPP variant index are only valid for SNVs: mirrored
    indels (T/TCA vs TCA/T) sort to the same key though they are distinct variants."""
    assert allele_cols, "assert_all_snv requires at least one allele column"
    non_snv = df.filter(
        pl.any_horizontal(pl.col(c).str.len_bytes() != 1 for c in allele_cols)
    )
    assert non_snv.height == 0, (
        f"{non_snv.height} non-SNV variant(s) in columns {allele_cols}; the unordered "
        f"allele key is invalid for indels. First rows:\n{non_snv.select(allele_cols).head()}"
    )
```

- [ ] **Step 4: Run the test + green + commit**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/ppp_database/test_allele_key_snv_guard.py -v` (PASS), then `pixi r invoke green`.

```bash
git add mecfs_bio/build_system/task/ppp_database/allele_key.py test_mecfs_bio/unit/build_system/task/ppp_database/test_allele_key_snv_guard.py
git commit -m "feat: add assert_all_snv guard and document unordered_allele_key as SNV-only"
```

---

### Task 2: SNV-assert the PPP HapMap3 variant index

**Files:**
- Modify: `mecfs_bio/build_system/task/ppp_database/construct_ppp_variant_index_task.py` (in `execute`, after `index` is collected, before `write_byte_stream_split_parquet`)
- Test: none (uses the Task-1 guard; the HapMap3 index is 100% SNV — verified — so this passes and drops nothing; the common-1kg mode is not built today).

**Interfaces:**
- Consumes: `assert_all_snv`, `GWASLAB_EFFECT_ALLELE_COL`, `GWASLAB_NON_EFFECT_ALLELE_COL`.

- [ ] **Step 1: Add the assertion**

Import `assert_all_snv` at the top and, immediately after `.collect()` produces `index`:

```python
        assert_all_snv(index, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL)
```

Add a one-line comment: the whole PPP database keys on the unordered {EA, NEA} set, which is invalid for indels, so the index must be SNV-only.

- [ ] **Step 2: Green + commit**

Run: `pixi r invoke green` (exercises PPP unit tests; HapMap3 index is SNV so the assert is satisfied).

```bash
git add mecfs_bio/build_system/task/ppp_database/construct_ppp_variant_index_task.py
git commit -m "feat: assert the PPP variant index is SNV-only (unordered key requires it)"
```

---

### Task 3: convert the baseline-LF annotation join path to exact (CHR, BP, A1, A2)

**Files:**
- Modify: `mecfs_bio/build_system/task/annotation_weights/build_baseline_lf_annotation_parquet_task.py` (`_dedup_one_chromosome`)
- Modify: `mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py` (`_JOIN_KEYS`, meta build ~122-130, annot join ~211-217)
- Modify: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py` (`_load_annotations`, `_ANNOT_KEY` usages, `_assert_annotation_keys_unique`)
- Test: none new; guarded by `test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py` (update it if it encodes the allele-key behavior).

**Interfaces:**
- The annotation matrix stays keyed/unique on exact `(CHR, BP, A1, A2)` and now retains mirrored-indel variants (restores the ~548 the unordered-key dedup dropped). Downstream joins match on the exact tuple; annotation `A1`→`NEA`, `A2`→`EA`.

- [ ] **Step 1: Producer — exact-key dedup**

In `_dedup_one_chromosome`, replace the unordered-key dedup with an exact one (remove the `unordered_allele_key` import and `_ALLELE_KEY_COL`):

```python
    lazy = pl.scan_parquet(member_path)
    schema = lazy.collect_schema()
    annot_cols = [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]
    key_cols = [_CHR_COL, _BP_COL, _A1_COL, _A2_COL]
    deduped = (
        lazy.with_columns([pl.col(c).cast(pl.Float32) for c in annot_cols])
        .unique(subset=[*key_cols, *annot_cols], keep="first")
        .collect()
    )
    conflicting_keys = deduped.height - deduped.n_unique(subset=key_cols)
    assert conflicting_keys == 0, (
        f"{member_path.name}: {conflicting_keys} (CHR, BP, A1, A2) group(s) carry "
        "differing annotations, violating the dedup assumption"
    )
    return deduped.sort(_BP_COL)
```

Update the module docstring: the result is unique on exact `(CHR, BP, A1, A2)`; mirrored indels are distinct rows (each allele orientation is its own variant). `_A1_COL`/`_A2_COL` already exist in this module.

- [ ] **Step 2: Ridge weights — exact join**

In `ridge_annotation_weights_task.py`: set `_JOIN_KEYS = [_CHR_COL, _BP_COL, _A1_COL, _A2_COL]`, delete `_ALLELE_KEY_COL` and the `unordered_allele_key` import. Change the meta build (~122-130) to keep A1/A2:

```python
            .select(_CHR_COL, _BP_COL, _A1_COL, _A2_COL, SNPVAR_COL)
            .unique(subset=[_CHR_COL, _BP_COL, _A1_COL, _A2_COL])
```

And the per-chromosome annot build (~211-217): remove the `unordered_allele_key` `with_columns`, keeping `annot_chrom = pl.scan_parquet(...).filter(...).collect()`, then `frame = annot_chrom.join(meta, on=_JOIN_KEYS, how="inner")`. Update the `_JOIN_KEYS` comment to say the join is exact on (CHR, BP, A1, A2), both sides reference-oriented (A1==REF).

- [ ] **Step 3: Explain contrast — map A1→NEA, A2→EA, join on `_KEY`**

In `polyfun_explain_contrast_task.py`: delete `_ALLELE_KEY_COL`, `_ANNOT_KEY`, and the `unordered_allele_key` import. In `_load_run_variants`, drop the `_ALLELE_KEY_COL` `with_columns` (the frame already has `_KEY`). In `_load_annotations`, replace the rename/key block with:

```python
    result = frame.rename(
        {
            _ANNOT_BP_COL: GWASLAB_POS_COL,
            _ANNOT_A1_COL: GWASLAB_NON_EFFECT_ALLELE_COL,  # annotation A1 == REF == gwas NEA
            _ANNOT_A2_COL: GWASLAB_EFFECT_ALLELE_COL,
        }
    )
    _assert_annotation_keys_unique(result, chrom, bp_min, bp_max)
    return result
```

Change the two `join(annot, on=_ANNOT_KEY, ...)` sites (~249, ~270) to `on=_KEY`, and in `_assert_annotation_keys_unique` change `annot.select(_ANNOT_KEY).n_unique()` to `annot.select(_KEY).n_unique()` (and its docstring reference).

- [ ] **Step 4: Green + update the contrast unit test if needed**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py -v`. If it constructs annotation/run frames assuming allele-key joins, update the fixtures to exact A1/A2 in reference orientation (A1==REF). Then `pixi r invoke green`.

Note: this invalidates the built annotation matrix + ridge weights (content changes: +~548 indels). They rebuild on next run; no code-cache issue.

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/build_system/task/annotation_weights/build_baseline_lf_annotation_parquet_task.py mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_contrast_task.py test_mecfs_bio/unit/build_system/task/polyfun_explain/test_polyfun_explain_contrast_task.py
git commit -m "feat: exact (CHR,BP,A1,A2) joins across the baseline-LF annotation path (indel-safe)"
```

---

### Task 4: SUSIE polyfun-prior join on exact (CHR, POS, EA, NEA)

**Files:**
- Modify: `mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py` (`align_data`, ~338-353; remove the now-unused `unordered_allele_key` import)
- Test: none new; guarded by `test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py` and the system test.

**Interfaces:**
- The prior join keys on `(CHR, POS, EA, NEA)`. `load_prior` already renames the prior's A1→NEA and A2→EA (and the prior is A1==REF), so the gwas and prior share orientation; the existing "prior does not cover N variants" check remains the coverage guard.

- [ ] **Step 1: Replace the prior join**

In `align_data`, replace the `if prior is not None:` block that builds `allele_key` on both sides with a direct exact join:

```python
    if prior is not None:
        n_before = len(joined)
        joined = joined.join(
            prior,
            on=[
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ],
            how="left",
        )
        missing = joined.filter(pl.col(_PRIOR_COL).is_null())
        if missing.height > 0:
            examples = missing.select(
                GWASLAB_CHROM_COL, GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL,
            ).head(5)
            raise ValueError(
                f"polyfun prior does not cover {missing.height} of {n_before} "
                f"(gwas intersect ld) variants; first missing:\n{examples}"
            )
        prior_out = joined[_PRIOR_COL].to_numpy()
    else:
        prior_out = np.ones(len(joined))
```

Remove the `from ...allele_key import unordered_allele_key` import (it is now unused in this module).

- [ ] **Step 2: Green + commit**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_susie_r_finemap_task.py -v` then `pixi r invoke green`. If a prior fixture in the unit test relied on orientation-agnostic matching, ensure its A1/A2 are reference-oriented.

```bash
git add mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py
git commit -m "feat: SUSIE polyfun-prior join on exact (CHR,POS,EA,NEA)"
```

---

### Task 5: palindrome-drop and chrom-range filter pipes

**Files:**
- Create: `mecfs_bio/build_system/task/pipes/drop_palindromes_pipe.py`
- Create: `mecfs_bio/build_system/task/pipes/chrom_range_filter_pipe.py`
- Test: `test_mecfs_bio/unit/build_system/task/pipes/test_drop_palindromes_pipe.py`

**Interfaces:**
- Produces: `DropPalindromesPipe(ea_col: str, nea_col: str)` — removes strand-ambiguous SNVs (A/T, T/A, C/G, G/C). `ChromRangeFilterPipe(chrom: int, start: int, end: int, chrom_col: str, pos_col: str)` — keeps rows with `chrom_col == chrom` and `start <= pos_col <= end`. Both are narwhals `DataProcessingPipe`s. These reproduce, as pipes, what `HarmonizeGWASWithReferenceViaAlleles` did inline, so Task 6 can drop that task without changing results.

- [ ] **Step 1: Write the failing pipe test**

```python
# test_mecfs_bio/unit/build_system/task/pipes/test_drop_palindromes_pipe.py
import narwhals
import polars as pl

from mecfs_bio.build_system.task.pipes.drop_palindromes_pipe import DropPalindromesPipe


def test_drops_only_strand_ambiguous_snvs():
    df = pl.DataFrame(
        {
            "EA": ["A", "C", "A", "G"],
            "NEA": ["T", "G", "G", "A"],  # rows 0 (A/T) and 1 (C/G) are palindromic
        }
    )
    out = (
        DropPalindromesPipe(ea_col="EA", nea_col="NEA")
        .process(narwhals.from_native(df).lazy())
        .collect()
        .to_native()
    )
    assert out.to_dicts() == [{"EA": "A", "NEA": "G"}, {"EA": "G", "NEA": "A"}]
```

- [ ] **Step 2: Run to confirm failure**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/pipes/test_drop_palindromes_pipe.py -v`
Expected: FAIL (ImportError).

- [ ] **Step 3: Implement both pipes**

```python
# mecfs_bio/build_system/task/pipes/drop_palindromes_pipe.py
"""Drop strand-ambiguous (palindromic) SNVs. Mirrors the palindrome handling that
HarmonizeGWASWithReferenceViaAlleles did inline, extracted so fine-mapping can apply it as a
gwas pipe once that task is no longer used. See is_palindromic_expr for the polars equivalent."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe

_PALINDROME_PAIRS = (("A", "T"), ("T", "A"), ("C", "G"), ("G", "C"))


@frozen
class DropPalindromesPipe(DataProcessingPipe):
    ea_col: str
    nea_col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        is_palindromic = None
        for ea, nea in _PALINDROME_PAIRS:
            cond = (narwhals.col(self.ea_col) == ea) & (narwhals.col(self.nea_col) == nea)
            is_palindromic = cond if is_palindromic is None else (is_palindromic | cond)
        return x.filter(~is_palindromic)
```

```python
# mecfs_bio/build_system/task/pipes/chrom_range_filter_pipe.py
"""Keep only rows in a single closed genomic interval. Mirrors the chrom-range filter that
HarmonizeGWASWithReferenceViaAlleles applied inline."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe


@frozen
class ChromRangeFilterPipe(DataProcessingPipe):
    chrom: int
    start: int
    end: int
    chrom_col: str
    pos_col: str

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(
            (narwhals.col(self.chrom_col) == self.chrom)
            & (narwhals.col(self.pos_col) >= self.start)
            & (narwhals.col(self.pos_col) <= self.end)
        )
```

- [ ] **Step 4: Run the test + green + commit**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/pipes/test_drop_palindromes_pipe.py -v` (PASS), then `pixi r invoke green`.

```bash
git add mecfs_bio/build_system/task/pipes/drop_palindromes_pipe.py mecfs_bio/build_system/task/pipes/chrom_range_filter_pipe.py test_mecfs_bio/unit/build_system/task/pipes/test_drop_palindromes_pipe.py
git commit -m "feat: add DropPalindromesPipe and ChromRangeFilterPipe"
```

---

### Task 6: remove `HarmonizeGWASWithReferenceViaAlleles` from the fine-mapping generators

**Files:**
- Modify: `mecfs_bio/asset_generator/fine_mapping_asset_generator.py`
- Modify: `mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py`
- Test: none new; guarded by the from-scratch system test `test_mecfs_bio/system/test_decode_me_polyfun_susie.py`.

**Interfaces:**
- Both generators stop building a `HarmonizeGWASWithReferenceViaAlleles` task. SUSIE receives `build_37_sumstats_task` directly, with a `gwas_data_pipe` = `CompositePipe([sumstats_pipe, UniquePipe(keep="none" on the 4 key cols), DropPalindromesPipe (only when palindrome_strategy == "drop"), ChromRangeFilterPipe (only when chrom_range is given)])`. `align_data`'s existing exact `(CHR,POS,EA,NEA)` join aligns the (reference-oriented) gwas to the (reference-oriented) LD panel — the removed task's allele-flip was a no-op for reference-oriented inputs (verified). This is indel-safe (exact join) and is what lets Task 8 restrict the harmonizer to SNVs.

- [ ] **Step 1: Build the replacement gwas pipe in `fine_mapping_asset_generator.py`**

Replace the `harmonized_sumstats_task = HarmonizeGWASWithReferenceViaAlleles.create(...)` block and every `gwas_data_task=harmonized_sumstats_task` with `gwas_data_task=build_37_sumstats_task` plus a shared `gwas_data_pipe` passed to each `SusieRFinemapTask.create(...)`:

```python
    gwas_pipe_steps: list[DataProcessingPipe] = [
        sumstats_pipe,
        UniquePipe(
            by=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL],
            keep="none",
            order_by=[GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL],
        ),
    ]
    if palindrome_strategy == "drop":
        gwas_pipe_steps.append(
            DropPalindromesPipe(ea_col=GWASLAB_EFFECT_ALLELE_COL, nea_col=GWASLAB_NON_EFFECT_ALLELE_COL)
        )
    if chrom_range is not None:
        gwas_pipe_steps.append(
            ChromRangeFilterPipe(
                chrom=chrom_range.chrom, start=chrom_range.start, end=chrom_range.end,
                chrom_col=GWASLAB_CHROM_COL, pos_col=GWASLAB_POS_COL,
            )
        )
    finemap_gwas_pipe = CompositePipe(gwas_pipe_steps)
```

Pass `gwas_data_task=build_37_sumstats_task, gwas_data_pipe=finemap_gwas_pipe` to each of the five `SusieRFinemapTask.create(...)` calls (base, strict, 1-cs, 2-cs — and keep `prior_info` as-is). Remove the `harmonized_sumstats_task` field from `BroadFineMapTaskGroup` and its return (or set it to `build_37_sumstats_task` if external code references it — grep first). Drop the now-unused `HarmonizeGWASWithReferenceViaAlleles` import; keep `ChromRange` (still a public param type).

- [ ] **Step 2: Same change in `polyfun_explain_fine_mapping_asset_generator.py`**

In `_build_shared_locus_inputs`, build the identical `finemap_gwas_pipe` and store the gwas task + pipe on `SharedFineMapInputs` (add a `gwas_data_pipe` field), then pass both to the `susie_uniform`/`susie_polyfun` `SusieRFinemapTask.create(...)` calls in `generate_polyfun_explain_group`. Remove the `HarmonizeGWASWithReferenceViaAlleles` usage/import.

- [ ] **Step 3: Green + commit (defer the system test to CI)**

Run: `pixi r invoke green`. Do NOT run the from-scratch system test locally — it runs in CI on the PR as the final result check (Task 6 + Task 4 are the result-affecting changes it guards). It exercises the **keep**-palindrome path (its lead variant is the palindrome chr1:173,855,298 {A,T}); the `drop` path is not separately guarded, which is acceptable. If CI later shows the lead variant moved, diff the harmonized variant set before/after for the chr1 locus (old `harmonize_task` output vs the new piped gwas∩ld set), reconcile, and — if the change is understood and intended — update the expected values in the test and **surface the proposed doc update to the user** (do not edit `docs/`).

```bash
git add mecfs_bio/asset_generator/fine_mapping_asset_generator.py mecfs_bio/asset_generator/polyfun_explain_fine_mapping_asset_generator.py
git commit -m "feat: fine-mapping generators consume harmonized sumstats directly (drop HarmonizeGWASWithReferenceViaAlleles)"
```

---

### Task 7: enable the SUSIE build-match assertion

**Files:**
- Modify: `mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py` (`SusieRFinemapTask.create`)
- Test: `test_mecfs_bio/unit/build_system/task/test_susie_build_match_guard.py`

**Interfaces:**
- Consumes: `FilteredGWASDataMeta.harmonization_info`, `HarmonizableReferenceTableMeta.harmonization_info` (both from Plan 1).
- `SusieRFinemapTask.create` asserts, at construction: the gwas meta is a `FilteredGWASDataMeta` with `harmonization_info` set, the ld-labels meta is a `HarmonizableReferenceTableMeta` with `harmonization_info` set, and their `build`s match. A mismatch or missing provenance raises `AssertionError`.

- [ ] **Step 1: Write the failing guard test**

```python
# test_mecfs_bio/unit/build_system/task/test_susie_build_match_guard.py
from pathlib import PurePath

import pytest

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
)


def _gwas(build: str) -> FakeTask:
    return FakeTask(meta=FilteredGWASDataMeta(
        id=AssetId("g"), trait="decode_me", project="me_cfs", sub_dir=PurePath("processed"),
        harmonization_info=HarmonizationInfo(build=build, ref_allele_col="NEA", pos_col="POS"),
    ))


def _ld_labels(build: str) -> FakeTask:
    return FakeTask(meta=HarmonizableReferenceTableMeta(
        group="ukbb_reference_ld", sub_group="chr1", sub_folder=PurePath("processed"),
        extension=".parquet", id=AssetId("ld"),
        harmonization_info=HarmonizationInfo(build=build, ref_allele_col="NEA", pos_col="POS"),
    ))


def test_build_mismatch_raises():
    with pytest.raises(AssertionError):
        SusieRFinemapTask.create(
            asset_id="s", gwas_data_task=_gwas("19"), ld_labels_task=_ld_labels("38"),
            ld_matrix_source=BroadInstituteFormatLDMatrix(FakeTask(meta=_ld_labels("38").meta)),
            effective_sample_size=1000,
        )
```

Match `SusieRFinemapTask.create`'s required args to its current signature (read it first); the point is a build-mismatched gwas/ld pair raises.

- [ ] **Step 2: Run to confirm failure**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/test_susie_build_match_guard.py -v`
Expected: FAIL (no assertion yet).

- [ ] **Step 3: Add the assertion in `create`**

Add a free helper near the module top and call it first in `create` (import `FilteredGWASDataMeta`, `HarmonizableReferenceTableMeta`):

```python
def assert_gwas_harmonized_to_ld(gwas_task: Task, ld_labels_task: Task) -> None:
    gwas_meta = gwas_task.meta
    ld_meta = ld_labels_task.meta
    assert isinstance(gwas_meta, FilteredGWASDataMeta) and gwas_meta.harmonization_info is not None, (
        "SUSIE gwas input must be a harmonized FilteredGWASDataMeta (run genome-reference harmonization)"
    )
    assert isinstance(ld_meta, HarmonizableReferenceTableMeta) and ld_meta.harmonization_info is not None, (
        "SUSIE ld_labels must be a HarmonizableReferenceTableMeta with harmonization_info"
    )
    assert gwas_meta.harmonization_info.build == ld_meta.harmonization_info.build, (
        f"gwas harmonized to build {gwas_meta.harmonization_info.build} but LD panel is build "
        f"{ld_meta.harmonization_info.build}"
    )
```

Call `assert_gwas_harmonized_to_ld(gwas_data_task, ld_labels_task)` at the start of `create`.

- [ ] **Step 4: Run the guard test + green + commit**

Run the guard test (PASS), then `pixi r invoke green`. Do not run the system test locally — in CI on the PR the real DecodeME gwas carries build "19" (Plan 1) and the LD labels build "19", so the guard passes there.

```bash
git add mecfs_bio/build_system/task/r_tasks/susie_r_finemap_task.py test_mecfs_bio/unit/build_system/task/test_susie_build_match_guard.py
git commit -m "feat: SusieRFinemapTask.create asserts gwas is harmonized to the LD panel's build"
```

---

### Task 8: SNV-assert the allele harmonizer + document remaining unordered-key uses

**Files:**
- Modify: `mecfs_bio/build_system/task/harmonize_gwas_with_reference_table_via_chrom_pos_alleles.py` (`execute`)
- Modify: `mecfs_bio/build_system/task/ppp_database/allele_key.py` (extend the docstring with the map of remaining uses — item 8)
- Test: none new; blast radius (LCV + the chr1_173 manual prep) is covered by green.

**Interfaces:**
- `HarmonizeGWASWithReferenceViaAlleles.execute` asserts both the gwas and the reference are SNV-only (its palindrome logic is SNV-only anyway). Other callers — `lcv_asset_generator.py`, the two LCV analyses, and `chr1_173_locus/harmonize_with_polyfun_reference_alleles.py` — must be SNV-only (LCV munges to HapMap3 SNPs, so it is); green surfaces any that are not.

- [ ] **Step 1: Add the SNV assertion**

In `execute`, after `gwas_data = _convert_ea_nea_to_str(gwas_data)` and `reference = _convert_ea_nea_to_str(reference)`, add (import `assert_all_snv`):

```python
        assert_all_snv(gwas_data, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL)
        assert_all_snv(reference, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL)
```

Update the class docstring: this task is valid only for SNVs (its palindrome detection is SNV-only); indel-containing data must use exact `(chrom, pos, ea, nea)` joins directly (as the fine-mapping path now does).

- [ ] **Step 2: Document the remaining unordered-key uses (item 8)**

Extend the `allele_key.py` module docstring with the current map:

```
Remaining unordered_allele_key uses, all SNV-safe by assertion or by construction:
  - construct_ppp_variant_index_task: SNV-asserted (HapMap3 index is SNV-only).
  - common_1kg_membership_task, build_slim_protein_parquet_task: the common-1kg PPP mode
    contains indels and is NOT built today; add filter-then-assert before enabling it.
The annotation path and the SUSIE prior/main joins no longer use this key (they join on the
exact (CHR, POS, EA, NEA) tuple).
```

- [ ] **Step 3: Green (blast-radius check) + commit**

Run: `pixi r invoke green`. If an LCV test or the chr1_173 prep trips `assert_all_snv`, that dataset carries indels — decide with the user whether to filter indels upstream for that path; do not silently weaken the assert.

```bash
git add mecfs_bio/build_system/task/harmonize_gwas_with_reference_table_via_chrom_pos_alleles.py mecfs_bio/build_system/task/ppp_database/allele_key.py
git commit -m "feat: restrict HarmonizeGWASWithReferenceViaAlleles to SNVs; document remaining unordered-key uses"
```

---

## Self-Review

**Spec coverage (Plan 2 = items 1-3, 5, 6, 8 + the SUSIE guard):**
- `assert_all_snv` helper + allele_key doc → Task 1. ✓
- PPP index SNV assert (item 2/3) → Task 2. ✓
- Annotation indel-aware rebuild, 3 coupled files (item 1) → Task 3. ✓
- SUSIE prior exact join (item 5) → Task 4. ✓
- Remove harmonizer from generators (item 6) → Task 6 (+ pipes in Task 5). ✓
- SUSIE build-match assertion (spec flow step 5) → Task 7. ✓
- Harmonizer SNV assert (item 3 / spec item "assert everything is SNV") → Task 8. ✓
- Document remaining uses (item 8) → Task 8. ✓

**Placeholder scan:** no "TBD"/"handle edge cases"; code steps show real code. "Read the current signature first" notes (SUSIE `create`, the contrast unit-test fixtures) are verification instructions.

**Type consistency:** `assert_all_snv(df, *cols)` used identically in Tasks 1/2/8. `DropPalindromesPipe(ea_col, nea_col)` / `ChromRangeFilterPipe(chrom, start, end, chrom_col, pos_col)` consistent in Tasks 5/6. Annotation `A1`→`NEA`, `A2`→`EA` mapping consistent in Task 3. `_KEY = [CHR, POS, EA, NEA]` reused (not `_ANNOT_KEY`) after Task 3.

**Trivial-test scan:** new tests only on `assert_all_snv` (guard), `DropPalindromesPipe` (filter), and the SUSIE build-match guard. Result-affecting changes lean on the existing from-scratch system + contrast tests. No wiring tests.

## Risks / notes for the executor
- **Result-affecting tasks are 4 and 6** (prior join, generator pipe change); the annotation rebuild (Task 3) affects the explainability tables/plots, not the finemap lead variant (the prior comes from the precomputed snpvar_meta, not the baseline-LF matrix). Verify Task 3 with the contrast unit test locally; verify Tasks 4 and 6 with the system test **in CI on the PR** (do not run it locally). If a result legitimately changes, update the expected values in the test and **surface the proposed doc update to the user** — never edit `docs/` (it is human-written).
- **Task 6 palindrome equivalence:** the `with_palindromes` loci pass `palindrome_strategy="keep"` + a `chrom_range`; the `without_palindromes` loci use `"drop"` + no range. The replacement pipe must reproduce both. The system test only covers the `keep` path (its lead variant is a palindrome), which is the primary case; the `drop` path is not separately guarded. If CI shows the lead variant moved, diff the pre/post harmonized variant set for the chr1 locus before touching expected values.
- **Task 8 blast radius:** the SNV assert also gates LCV and the chr1_173 manual prep — both expected SNV-only, but confirm via green rather than assuming.
- Plan 1 must be merged first (Task 7 needs `harmonization_info`; Task 6 assumes the LD labels are a `RenameColsTask` output carrying it).
