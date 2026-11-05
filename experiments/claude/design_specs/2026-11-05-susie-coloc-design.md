# SuSiE-coloc for trait × UKB-PPP colocalization: design

Date: 2026-11-05. Status: draft for review.

## Intent

The merged coloc.abf pipeline (TraitProteinColocTask, PR #1327) assumes at most one causal variant
per trait and per protein in each locus window. When a locus carries several signals, that
assumption can force H3 (distinct causal variants) even when one of the signals is shared. In the
four v1 runs, 93 of 98 tested pairs ended at PP_H3 ≈ 1. A screen of those results
(experiments/claude/susie_coloc_probes/abf_h3_multisignal_screen.py) found pairs where the protein's
lead variant is itself strongly associated with the trait but is not the trait lead:

| Trait | Protein | Trait −log10P at its lead | Trait −log10P at the protein's lead | Distance between leads | abf PP4 |
|---|---|---|---|---|---|
| LDL | APOE | 742 | 160 | 2 kb | 0 |
| LDL | ABO | 49 | 41 | 0.7 kb | 0.027 |
| IBD | IL12B | 39 | 31 | 56 kb | 5e-6 |

SuSiE-coloc (Wallace 2021) fine-maps each side with SuSiE and runs coloc on every pair of signals,
which is the standard fix for this. This project adds it as a second colocalization task that runs
alongside coloc.abf. It does not replace coloc.abf.

## Scale

Up to about 10 traits × 50 proteins. That means roughly 30–50 SuSiE-coloc pairs per trait, spread
over about 2–15 tested loci per trait. The v1 abf runs tested 2, 9, 10 and 14 distinct loci.

## Goals

- A new Task, separate from TraitProteinColocTask, that writes SuSiE-coloc results for one trait
  against the UKB-PPP proteins selected at its loci.
- Reuse the abf pipeline's locus finding, candidate selection and budget, dense protein fetching,
  protein orientation and exact alignment.
- Reuse the existing SuSiE fine-mapping machinery (susieR through rpy2, the LD-repair and
  diagnostic helpers in susie_r_finemap_task.py, and the Broad UKBB LD URL helpers).
- Wire it into the four traits that already run coloc.abf: Willer LDL, DecodeME, deCODE
  seropositive RA and Liu et al. 2023 IBD.

## Non-goals

- A pure-Python SuSiE. The probes below show SuSiE is not the bottleneck. Keep susieR.
- Any on-disk cache of LD matrices, whether in the asset store or elsewhere.
- Replacing or changing coloc.abf or its outputs.
- SuSiE's refinement procedure, and susieR 0.16's R_finite / R_mismatch model. Only 0.14.2 is on
  conda-forge.
- Display tables and locus plots for the new results. Revisit after v1 results are in.
- Falling back to coloc.abf when SuSiE finds no credible set.

## Anti-goals (failure even if it technically works)

- Disk use that grows with the number of loci or traits.
- Losing PCSK9 × LDL colocalization.

## Probe results that shaped the design (2026-11-05)

- Fetching one dense UKB-PPP protein from Synapse took 17–22 s, the same as in the abf runs.
- A Broad UKBB LD matrix covers 3 Mb in 1 Mb steps. The median file is 1.08 GB, the 90th
  percentile 1.56 GB and the largest 4.7 GB; all files together are 3 TB. Downloading a
  median-size file took 85 s (about 12 MB/s), and scipy.sparse.load_npz took 17 s. The npz members
  are deflate-compressed COO arrays, so S3 range reads cannot fetch a subset.
- Rerunning the base DecodeME chr1 fine-mapping task (2,009 variants) took 55 s: about 20 s of
  start-up, 20 s of loading and aligning, 4 s of PSD repair, 10 s for the estimate_s_rss and
  kriging_rss diagnostics, and 4 s for susie_rss.
- susie_rss on slices of a real LD matrix with simulated z-scores and L = 10:

  | p | numpy-to-R copy | susie_rss |
  |---|---|---|
  | 1,000 | 0.5 s | 1.1 s |
  | 2,000 | 2.0 s | 1.5 s |
  | 4,000 | 7.0 s | 4.0 s |
  | 6,000 | 14.6 s | 9.3 s |

- Estimated cost per trait: about 15–20 min of protein downloads, 15–25 min of LD downloads, and
  under 20 min of SuSiE including the copies.

Scripts are in experiments/claude/susie_coloc_probes/ and logs in experiments/claude/logs/.

## Decided by the user

1. LD matrices are downloaded on the fly, used, and discarded. Nothing is cached.
2. A pair where SuSiE finds no credible set on the trait side or the protein side gets one row
   with status untestable and no posteriors. There is no abf fallback.
3. Results have one row per (trait credible set, protein credible set) combination, so a pair can
   have several rows. Filtering happens downstream.
4. Success: PCSK9 × LDL must colocalize. Whether APOE, ABO and IL12B find a shared signal is
   reported, not required.
5. Use the out-of-sample-LD diagnostics and filtering from the existing fine-mapping task, on both
   sides.
6. Before running diagnostics, drop variants whose N is below 80% of the largest N in that side's
   window.
7. A separate task, alongside abf.

## Design

### Flow for one trait

1. Find loci and select (locus, protein) pairs exactly as TraitProteinColocTask does: find the
   trait loci, build candidates from the hit table, then select within the budget. Statuses for
   excluded MHC, not selected by budget, and no candidates carry over unchanged.
2. Fetch, orient and align every selected protein's windows as the abf task does. These are small
   and are held for the whole run.
3. Loop over tested loci, one at a time:
   1. Choose the Broad LD block for the locus window, and download its labels (.gz) and matrix
      (.npz) into the task's scratch directory.
   2. Load the matrix, keep only the rows and columns of the variants used at this locus,
      densify, then delete both files and free the full sparse matrix before any SuSiE work.
   3. Fine-map the trait once on this locus's variants (see "Per-side SuSiE").
   4. For each protein paired with the locus, fine-map the protein, then run coloc on every
      combination of trait credible set and protein credible set (see "Combining").
   5. Discard the LD submatrix.

At most one LD file is on disk at a time. Peak memory is about 3 GB while a file loads (215M
non-zeros: int32 row and column indices plus float32 data).

### Mapping our variants to the hg19 LD matrix

The Broad matrices are on GRCh37. Our trait and protein windows are on GRCh38. UKB-PPP's regenie ID
column encodes the GRCh37 position. The PPP variant index already parses this as POS_HG19, with no
liftover. Every variant passed to colocalization is in both trait and protein, so every one has a
PPP hg19 position.

- Extend the dense window fetch to keep the hg19 position parsed from ID. Carry it through
  orientation and alignment.
- Match aligned variants to the Broad labels on (hg19 position, unordered allele pair).
- Orient every variant's z-scores on both sides so the effect allele is the Broad label's allele2.
  The sign convention does not matter as long as every variant uses the same one: flipping every
  variant leaves R unchanged. What must hold is that each variant's z sign matches the coding of
  its row and column.
- Variants without a Broad match are dropped and counted.

### Choosing the LD block

Blocks are 3 Mb wide and start every 1 Mb, so any window up to 2 Mb fits inside some block. Merged
windows can be wider; the LDL APOE window is 2.14 Mb. Choose the block that contains the most of
the locus's aligned variants (ties go to the most central). Variants outside the block are dropped
and counted, and the row records how many. The existing get_optimal_ukbb_ld_interval chooses by a
single position, so it needs a variant-count variant.

### Per-side SuSiE

The same steps run on the trait side (once per locus) and on each protein side (once per pair):

1. Drop variants with N < 0.8 × the side's maximum N in this window. If a side has no per-variant
   N (it has a single sample size), this step drops nothing.
2. Repair the LD submatrix to be positive semi-definite (make_psd_corr).
3. Estimate λ with estimate_s_rss.
4. Run kriging_rss and drop variants with logLR > 2 and |z| > 2.
5. Run a second round, as Zou et al. 2022 recommend: re-estimate λ and rerun kriging on the
   filtered variants, dropping again with the same thresholds.
6. Regularize to (1 − λ)R + λI, using the λ from the last round.
7. Run susie_rss(z, R, n, L = 10) with the residual variance fixed (the default).
8. For n, use the median N of the variants that remain.

The trait and protein can end with different variant sets after filtering. Combining uses their
intersection.

Implementation detail for the builder: the numpy-to-R copy of R costs as much as susie_rss at large
p. Use a cheaper transfer if one exists, but this is not required.

### Combining (coloc.susie logic, in Python)

For every credible set on each side, take that SuSiE effect's per-variant log Bayes factors
(lbf_variable). Restrict both to the variants they share. Compute PP0–PP4 with the existing
combine_abf in coloc_abf.py, using the same priors as the abf runs (p1 = p2 = 1e-4, p12 = 1e-5).
coloc.susie's own default is p12 = 5e-6; my call is to match our abf runs so the two are
comparable. Only effects that form a credible set are used, as in coloc.susie.

### Statuses

- tested: one row per credible-set combination.
- untestable: SuSiE found no credible set on the trait side or the protein side. One row, with
  empty credible-set and PP columns.
- failed with a reason: too few variants after filtering (reuse min_aligned_variants), SuSiE did
  not converge, or a download error.
- excluded MHC, not selected by budget, no candidates: as in abf.

### Output (raw_results.parquet)

All abf identification and selection columns, plus:

- trait and protein credible-set indices;
- for each credible set: lead variant (position and alleles), size, minimum purity, and the
  SuSiE effect's log Bayes factor;
- PP_H0 to PP_H4 and the priors;
- per side: λ from each round, variants dropped for low N, variants dropped by kriging in each
  round, variants without a Broad label match, variants outside the LD block, variants passed to
  SuSiE, and n;
- the Broad block used.

Integers are stored as Int32, as in abf.

### Code organization

- New task in mecfs_bio/build_system/task/protein_coloc/.
- Move steps 1–2 out of TraitProteinColocTask into shared free functions, so both tasks call the
  same code rather than duplicating it.
- LD access goes behind a small Protocol (a source that takes a locus's variants and returns an
  aligned LD submatrix, plus labels), with a Broad UKBB implementation. This mirrors
  DenseWindowFetcher and lets unit tests inject a fake without mocking. Move the Broad
  block-selection and URL helpers it needs from asset_generator/ukbb_broad_ld_matrix_generator.py
  into build_system, because build_system cannot import from asset_generator.
- Download with wf.download_from_url into scratch, and delete immediately after loading.
- Extract the SuSiE preprocessing (PSD repair, λ, kriging filter, regularization) from
  susie_r_finemap_task.py into shared functions, rather than copying them.
- Wire one new task per trait next to the existing abf task in ppp_coloc_asset_generator.py.

### Testing

- A task-level unit test on synthetic data with a fake LD source and a fake window fetcher. It
  should cover a shared signal giving high PP4, a no-credible-set side giving untestable, and the
  80% N filter. R-backed tests go in their own module.
- Validation runs on the four traits:
  - pass/fail: PCSK9 × LDL has a credible-set combination with PP4 > 0.8;
  - also checked: CELSR2, MST1 and CCL19;
  - reported, no threshold: APOE, ABO and IL12B.

## Smaller calls made while writing this (Claude's, not the user's)

- The trait must have an N column. Assert this at the start, because SuSiE needs n.
- Case-control traits use z = beta / se and their effective N. No special handling.
- Credible sets use susieR's defaults: 95% coverage and minimum purity 0.5.
- The pair budget and scoring use the same SelectionConfig defaults as the abf task, and can be
  overridden per trait.
- LD downloads are retried with call_with_retries.

## Left to the builder

- The exact Protocol shapes and function boundaries.
- How the hg19 position is threaded through orient_and_align.
- Faster numpy-to-R transfer, if one is easy.
- Whether to save per-locus SuSiE PIPs or aligned windows for later plotting. They are cheap, and
  the abf task saves aligned windows.
