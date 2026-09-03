# Polyfun-explain annotation callouts — design

Date: 2026-08-30
Branch: polyfun-precomputed-prior

## Goal

Replace the dropped annotation-family panels with a small number of on-figure
callouts that name, for the variants the polyfun prior most clearly boosted,
which annotation families drove that boost. One callout per polyfun credible
set, at most, e.g.

    173845678:A:T (conserved ++, coding +)

read as `pos:nea:ea` followed by the key annotation families with a `+`/`++`
strength marker.

The callout answers "which variants did the functional prior single out here,
and on what functional grounds" without a genome-wide family track that was too
jagged to read.

## Where the work lives

All selection and text is computed in `PolyfunExplainContrastTask`
(deterministic, table-in/table-out, unit-testable). It emits a new
`callouts.parquet`. `PolyfunExplainPlotTask` reads that table and only places
the labels — it derives no statistics of its own. The current single-focal
`selection.json` becomes the degenerate case of this table; keep it for now (the
plot no longer reads it) or fold it in, but do not regress its tests.

The contrast task already has everything needed: both runs' variants and CS
tables, the full per-variant PIP files (`PIP_FILENAME`), gamma (`GAMMA_RAW_COL`),
the annotation matrix, `abar_c` (uniform-PIP-weighted annotation means), and the
per-family contrast `C_family(i) = family_scaled(i) − Σ_c γ_c·abar_c`.

## Callout selection (which variants)

Consider each polyfun credible set independently. Within a CS, take the variant
with the highest polyfun PIP as the candidate. Flag it iff BOTH gates pass:

1. **Dominance.** Its polyfun PIP is at least 0.05 (5 pp) above the
   next-highest polyfun PIP in the same CS. For a single-member CS the "next"
   PIP is 0. (This gate already guarantees at most one flagged variant per CS:
   only the max can be ≥5 pp above the second-max. Exact ties → no flag, read as
   "unresolved".)

2. **Prior effect.** Its polyfun PIP is at least 0.10 (10 pp) above the same
   variant's PIP in the **uniform** run, read from the uniform run's full
   per-variant PIP file (`PIP_FILENAME`), NOT uniform CS membership. A variant
   absent from the uniform run → uniform PIP 0.

No absolute PIP floor: selection is purely about the prior-induced change, per
the feature's goal.

Variant identity across runs and against the annotation table uses the existing
allele-aware key (`_KEY` = CHR, POS, EA, NEA; annotations joined on the
unordered-allele key), so no new matching logic.

## Key families per callout (which annotations, and strength)

Computed per flagged variant (NOT one global top-3), positive side only:

1. For every family, the "diff" is the existing per-family contrast at that
   variant, `diff_f(i) = family_scaled_f(i) − Σ_{c∈f} γ_c·abar_c` (already
   computed in `per_family`).

2. Normalize by the spread of the family's scaled value over the background:
   `sd_f` = the **uniform-PIP-weighted** standard deviation of `family_scaled_f`
   over the **same** variant set and weights used for `abar` (uniform-run
   variants, uniform PIP weights). Using the uniform-PIP background means the
   flagged variant's own weight in mean/SD is its tiny uniform PIP, so no
   leave-one-out correction is needed.

3. Keep only families with `diff_f(i) > 0` (elevated at this variant — a
   negative contrast is a family arguing against prioritization and is not part
   of "why it was boosted").

4. Bucket by `z_f = diff_f(i) / sd_f`:
   - `z_f ≤ 1`: drop (not distinguishing).
   - `1 < z_f ≤ 2`: mark `+`.
   - `z_f > 2`: mark `++`.

5. Of the surviving families, keep the top 3 by `z_f`, ordered by `z_f`
   descending, and render as `family_short_label marker`, comma-separated,
   using `FAMILY_SHORT_LABELS`.

Edge cases:
- `sd_f ≈ 0` (family constant across the locus): skip that family (no divide).
- No family survives: render just `pos:nea:ea` with no parentheses.

## `callouts.parquet` schema

One row per flagged variant:

| column      | meaning                                             |
|-------------|-----------------------------------------------------|
| CHR/POS/EA/NEA | the existing `_KEY` columns                      |
| `cs`        | polyfun credible-set number                         |
| `pip_pf`    | polyfun PIP (the anchor y for the label)            |
| `pip_u`     | uniform PIP of the same variant                     |
| `label`     | the full rendered string, e.g. `173845678:A:T (conserved ++, coding +)` |

Deterministic and independent of rendering, so it can be asserted directly in a
Task-level test.

## Plot rendering

`PolyfunExplainPlotTask` reads `callouts.parquet` and places one label per row
on the **PIP (polyfun)** panel, anchored at `(POS, pip_pf)`, using `textalloc`
(`ta.allocate`) as in
`mecfs_bio/build_system/task/magma/magma_plot_brain_atlas_result_with_stepwise_labels.py`
— pass the panel's stem points as the scatter to avoid, with leader lines to the
anchor.

- Give the polyfun PIP panel headroom so labels have somewhere to land: raise its
  ylim above the shared PIP top, or widen its GridSpec row via `height_ratios`.
  Keep the two PIP panels' *data* y-scale shared for comparability even if the
  polyfun panel's drawn ylim is taller.
- Verify `ta.allocate` determinism (seed / stable placement) so repeated builds
  don't churn the SVG. L=1 has ≤1 label (trivial); L=10 is the stress case.
- Empty `callouts.parquet` → render exactly today's figure (no labels).

## Testing

Task-level tests on the contrast task, on the existing synthetic explain inputs
(`build_synthetic_explain_inputs`), rigged so:
- a variant clears both gates and yields ≥1 family past 2 SD (`++`) and one
  between 1–2 SD (`+`) → assert its `label`;
- a CS with a near-tie top-2 (dominance fails) → no row;
- a variant high in both runs (prior-effect gate fails) → no row;
- a flagged variant with no family past 1 SD → label has no parentheses.

Assert on the `callouts.parquet` contents (values, not message text). The plot
task keeps only its existing "writes png+svg" smoke test.

## Decisions (settled)

1. No absolute PIP floor — change-based only.
2. Key families computed per callout, not one global top-3.
3. Positive-side families only.
4. Callouts computed in the contrast task; plot reads `callouts.parquet`.
