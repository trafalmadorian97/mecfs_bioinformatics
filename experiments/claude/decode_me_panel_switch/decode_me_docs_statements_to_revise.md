# DecodeME docs: statements invalidated by the regenerated figures

This lists statements in `docs/` that the regenerated DecodeME figures no longer support, for
rewriting by hand. The comparison is against the figures that were published before the
regeneration. Those figures predate both the gnomAD panel switch (PR #1233) and the 2026-10-07
baseline (`decode_me_panel_switch_results.md`). Some small numeric changes, such as the LDSC
intercept moving from 0.9421 to 0.9423, therefore come from earlier changes, not from the panel.

## How the figures were refreshed

| Step | Script | Result |
|---|---|---|
| Export every DecodeME figure task that depends on the gnomAD panel, and hash each file against the manifest | `find_stale_figures.py` | 92 tasks; 81 files changed in bytes, 11 unchanged |
| Separate content changes from incidental byte changes (plotly ids, PNG metadata) | `classify_figure_changes.py` | 80 changed in content; 1 changed in bytes only (curated gene-set MAGMA bar plot, left as is) |
| Delete the 80 manifest entries and local files, then regenerate them into `docs/_figs` with `generate_new_figures` | `refresh_decode_me_figures.py` | Done; the manifest now lacks these 80 entries |
| Summarize old vs new content | `describe_figure_changes.py` | Source of everything below |

Out of scope:
- Univariate MiXeR figures: not regenerated, because refitting is slow. Their published versions
  are kept.
- `decode_me_minus_*` subtraction figures, multi-trait figures and other studies: deferred.

Still to do:
- Run `push_figures` (or `invoke publish-figures`) to re-add the 80 entries with their new hashes
  and upload the blobs.
- Commit the manifest.

## Statements that are now wrong

### SuSiE (uniform prior), chr20 47M–48.2M

`docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE/f_Chr20_47M_48M_Locus.md`

- **Line 18:** "the $L=10$ run returned three credible sets, while the $L=2$ run returned two …
  One of these credible sets is highly concentrated … while the others are diffuse."
  - The $L=10$ run now returns **two** credible sets, the same as $L=2$:
    - L2 is unchanged: 20:47743125:C:A, PIP 1.0.
    - L3 is diffuse: 86 variants, top 20:47731228:A:C with PIP 0.107.
  - The old diffuse L1 (119 variants, spanning roughly 47.51–47.91 Mb) is gone. Its variants
    are mostly absorbed into L3 (79 of 119 shared).
  - "The others are diffuse" should become "the other is diffuse".
- **Line 21:** "the variants found by the $L=1$ and strict $L=10$ runs are identical, while the
  other runs differ."
  - The first half still holds: $L=1$ and strict $L=10$ share the same 58 variants.
  - The new UpSet plot also shows that **$L=10$ and $L=2$ are identical**: both have 87 variants.
    The plot now has two groups instead of four runs that each differ:
    - 42 variants are shared by all four runs;
    - 45 are only in $L=10$ and $L=2$;
    - 16 are only in $L=1$ and strict $L=10$.
- **Line 31:** "the highly concentrated nature of one of the credible sets and the diffuse nature
  of the others."
  - There is now one concentrated set and one diffuse set.

### SuSiE (uniform prior), chr17 50M–51M

`docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE/e_Chr17_50M_51M_Locus.md`

- **Line 23:** "The distribution of probability mass across these variants was slightly different
  across the runs, as illustrated by the following UpSet plot, which is restricted to variants of
  $PIP>0.01$."
  - The regenerated PIP > 0.01 UpSet plot shows the **same 27 variants in all four runs**.
  - In the published plot, $L=1$ lacked one of the 26 variants that the other runs had.
  - The plot no longer shows a difference between runs.

### PolyFun SuSiE (external prior), chr17 50M–51M

`docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE-PolyFun_(External_Prior)/e_Polyfun_Chr17_50M_51M_Locus.md`

- **Line 52:** "focusing especially on the three evolutionarily conserved variants
  **17:50291040:C:T**, **17:50299079:G:A**, and **17:50296459:G:T**."
  - The first two are still the top variants:

    | Variant | Old PIP | New PIP |
    |---|---|---|
    | 17:50291040:C:T | 0.229 | 0.226 |
    | 17:50299079:G:A | 0.185 | 0.183 |
    | 17:50296459:G:T | 0.101 | 0.099 |

  - The third variant is still third, but it has dropped just below 0.1.
- **Line 55:** "The next table shows the functional annotations for these top variants."
  - The per-variant annotation table now has **only two variant columns**: 17:50291040:C:T and
    17:50299079:G:A.
  - 17:50296459:G:T has dropped out, since its PIP fell below the table's cut-off. Text that
    names three variants in that table no longer matches it.

### Harmonization reference (not a figure, but no longer true after the switch)

These lines say palindromes and ambiguous indels were oriented using the Thousand Genomes
Project. DecodeME harmonization now uses gnomAD v2.1.1 genomes, the `nfe_nwe` (North-Western
European) group, with a palindrome allele-frequency distance cap of 0.1.

- `docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE/a_Chr1_173M_174M_Locus.md`, line 29
- `docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE-PolyFun_(External_Prior)/a_Polyfun_Chr1_173M_174M_Locus.md`, line 26
- `docs/Analysis/ME_CFS/DecodeME/Fine_Mapping/SUSIE-PolyFun_(Internal_Prior)/a_Polyfun_Chr1_173M_174M_Locus.md`, line 32

## Changed figures whose surrounding text is still correct

I checked these against the new figures and found no statement that needs rewriting. Some of
the numbers in the figures have moved.

| Page | What changed in the figures | Claims still true |
|---|---|---|
| `f_LDSC-DecodeME.md` | h2_liab 0.0814 → 0.08143; intercept 0.9421 → 0.9423; diagnostic plot bins identical | "0.0814" (still rounds to it); intercept < 1; bins below χ² = 1 |
| `h_MAGMA-HBA-DecodeME.md` | Independent-cluster p-values shift in the 4th significant figure | Same 3 independent clusters, in the same order |
| `i_MAGMA-H_DecodeME.md` | Newer plots draw only genes with −log10 p > 2; DNAH10 rises in the iPSC astrocyte map | Bonferroni-significant gene sets identical in all 6 maps; BTN1A1 still significant |
| `j_S-LDSC_DecodeME_Analysis.md` | Plot HTML only (plotting code); data points identical | All |
| `e_Cis-PQTL-MR-DecodeME.md` (region plots) | A few extra points (newly harmonized variants); lead variants and LD pattern unchanged | Non-colocalization reading of both plots |
| SuSiE chr1 (`a_…`) | Credible set 76 → 75 variants | Single diffuse set, identical across all 4 runs |
| SuSiE chr6 97M (`c_…`) | New top variant 6:98432302:C:CA (an indel, PIP 0.043), previously absent; set still 84 variants | Single set per run, "largely but not entirely overlapping" |
| SuSiE chr15 (`d_…`) | 149 → 147 variants; 15:55158922 PIP 0.194 → 0.195 | Identical across runs; concentrated on 55158922 |
| SuSiE chr17 (`e_…`) | 62 → 63 variants | Single set identical across runs (line 23 aside, see above) |
| SuSiE chr20 strict L=10 | 56 → 58 variants; same top variants | Single diffuse set |
| PolyFun external chr1 | 28 → 26 variants; top variants and PIPs unchanged to 2 decimals | All |
| PolyFun external chr15 | 149 → 147 | All |
| PolyFun external chr17 | Credible-set counts and UpSet plots: L=1 still lacks one variant | "L=1 run differs slightly" (lines 52 and 55 aside, see above) |
| PolyFun external chr20 | PolyFun sets unchanged (2 sets: 1 + 35 variants); the uniform column in the explain plot and table now shows 2 sets, not 3 | 20:47743125:C:A at PIP 1.0 in L=10 and L=2; 20:47663382:A:G has substantial PIP in every run |
| PolyFun internal chr1 | Max lift_ext 29.07 → 29.24; max lift_int 8.49 → 8.54; different top variants under each prior | All |
| PolyFun internal chr15 | 150 → 148 | All |
| PolyFun internal chr17 | The 50% PIP UpSet plot used to be identical across runs; now L=1 differs by one variant; 17:50291040:C:T uniform PIP 0.012 | "very similar, L=1 differs slightly"; "about 1 percent" |
| PolyFun internal chr20 | L=1 set 52 → 55 variants; 20:47526493:T:G still the top internal-prior variant (0.161) | All |

## Reproducing

```
pixi r python -m experiments.claude.decode_me_panel_switch.find_stale_figures
pixi r python -m experiments.claude.decode_me_panel_switch.classify_figure_changes
pixi r python -m experiments.claude.decode_me_panel_switch.refresh_decode_me_figures
pixi r python -m experiments.claude.decode_me_panel_switch.describe_figure_changes
```

`describe_figure_changes` downloads the published versions from the figure release into
`experiments/claude/decode_me_panel_switch/published_figs/` (not committed). It reads the
removed entries from the manifest at git HEAD, so it must run before the manifest change is
committed.
