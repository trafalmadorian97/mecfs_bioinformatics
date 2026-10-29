# Large-scale trait × proteomics colocalization — design

Date: 2026-10-29. Status: design, not yet planned or built.

## Intent

We already run large-scale genetic correlation of one trait against every UKB-PPP protein
(`mecfs_bio/build_system/task/ppp_ldsc/`). Genetic correlation is genome-wide and
underpowered for sparse signals; colocalization asks a sharper, locus-level question: at each
of the trait's GWAS hits, does any protein's pQTL share the same causal variant? A positive
answer nominates the protein as a candidate mediator of the trait association, which is the
kind of result that feeds directly into the DecodeME docs.

The pipeline is meant to be run on roughly 5–10 trait GWAS (more than once, but not many
times), so it fetches dense protein data on demand and discards it, rather than maintaining a
permanent dense store.

## Goals

1. Given one trait GWAS, produce colocalization results between each of its genome-wide
   significant loci and every UKB-PPP protein that has a significant hit in the same locus.
2. Outputs per trait run:
   - (a) a **raw results table** with the gory details, one row per (trait locus, protein) pair;
   - (b) a **display table** of high-probability colocalizations only (protein, location,
     H4), embedded in the docs with the `data_table` macro (`docs_macros/macros.py`);
   - (c) for each high-probability colocalization, a **stacked locus plot**: trait
     Manhattan, protein Manhattan, gene track, on a shared genomic x-axis.
3. **Source-agnostic.** Nothing in the colocalization core knows about PPP. The proteomics
   dataset plugs in through a small source interface so that the Western et al. 2024 CSF
   pQTLs (branch `csf-database`, unmerged) or another proteomics GWAS can be added later.
4. **Method-extensible.** v1 runs coloc.abf. The colocalization method sits behind its own
   interface so SuSiE-coloc can be added later as a second implementation without
   restructuring the pipeline.

## Non-goals (v1)

- CSF source (follow-up; see "Follow-ups").
- SuSiE-coloc and anything requiring LD (follow-up). Consequently, no LD download and no LD
  r² colouring in plots.
- A permanent dense PPP store on `/mnt/d`, a region-major store, or the VIPRS int8 LD
  matrices. Revisit only if usage grows well beyond ~10 traits.
- Mendelian randomization, PWAS, or any causal-direction inference.

## Anti-goals (failure even if it "works")

- **Silently joining distinct variants.** An exact join on (CHR, POS, EA, NEA) is *not*
  fail-safe if either side has a mis-oriented ambiguous indel: a deletion TCA→T stored flipped
  as (NEA=T, EA=TCA) collides with a genuine insertion T→TCA. Orientation must be established,
  not assumed (see "Allele orientation").
- **Unordered allele keys.** They conflate mirrored pairs (T/TCA vs TCA/T are distinct
  variants); this repo is actively retiring them. Do not use `unordered_allele_key` here.
- **Silent build mismatch.** Every join happens in GRCh38 with build recorded in metadata
  and asserted.
- **Overstating trans results.** Trans-pQTL hotspots (ABO, APOE, SH2B3, HLA, complement)
  colocalize with many proteins via one pleiotropic variant; every result must carry a
  cis/trans tag so readers can discount trans hits.

## Constraints

- 16 GB laptop runner (intentional regression guard). Memory per step must be bounded by a
  locus window or a chromosome, never a whole genome-wide protein file.
- Each PPP protein raw upload is a ~545 MB `.tar` on Synapse containing one gzipped regenie
  file per chromosome (~10–35 MB each). Sustained download ≈ 40–60 MB/s.
- Repo conventions apply (polars, attrs frozen, Task `.create()` derives meta from deps,
  inject dependencies instead of mocking, `execute()` reads as a high-level spec, etc.).

## Decisions made with the user

### Genome build: GRCh38 throughout

The analysis build is GRCh38. PPP `GENPOS` is hg38 (the regenie `ID` string carries hg19).
Western CSF is GRCh38. DecodeME's raw form is GRCh38. Traits only available in GRCh37 must be
lifted to 38 and then harmonized before use. (When SuSiE-coloc arrives, the GRCh37 Broad UKBB
LD matrices are used unchanged; only their variant labels are lifted 37→38, dropping variants
that fail liftover, change chromosome, or whose REF changes. That is a follow-up problem.)

### Trait side: must be genome-reference harmonized in GRCh38

The trait input must have been produced by `GenomeReferenceHarmonizationTask`, i.e. its meta
(`FilteredGWASDataMeta`) carries `harmonization_info` with `build == "38"`. Assert this in the
pipeline Task's `.create()` so a mismatch fails at graph construction. Trusted-path trait
harmonization (ambiguous indels keep source orientation) is accepted: trust there is earned
by full reference consistency on checkable variants.

Note: no GRCh38-harmonized DecodeME asset exists yet on main
(`decode_me/processed_gwas_data/build_38/` only assigns rsIDs). Wiring one
(`GenomeReferenceHarmonizationTask` with `UCSC_HG38_INDEXED_FASTA` and
`GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES`) is a prerequisite asset for the first real run.

### Protein side: run each locus window through the untrusted harmonization path, on the fly

Probe results (`experiments/claude/coloc_allele_orientation/probe.py`, logs
`experiments/claude/logs/coloc_allele_orientation_probe_*.log`):

- **PPP** (BTN1A1, RABGAP1L): vs hg38 FASTA at GENPOS, SNVs are 99.74% NEA(ALLELE0)==REF,
  0.16% EA==REF, 0.10% neither. 99.4% of the non-conforming SNVs have ALLELE0 == hg19 REF at
  the `ID` position, so ALLELE0 is the hg19 VCF REF and hg38 mismatches come from build drift.
- **CSF** (3 aptamers): positions perfect (0% neither), but effect_allele is the in-sample
  *minor* allele (15.5% of SNVs have EA==REF).
- Indels: ~92–93% have both alleles matching the reference prefix ("ambiguous indels");
  orientation of these cannot be established from the genome alone.

Decision: for every protein window, call the existing pure function
`resolve_chromosome(rows, ChromosomeContext(trusted=False, ...), load_panel)` from
`mecfs_bio/build_system/task/genome_reference_harmonization/resolve_chromosome.py`. This
orients SNVs by genome lookup, resolves ambiguous indels and palindromic SNVs by comparing
source EAF against panel AF, and drops anything indecisive. The panel is
`GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES` (European/NFE AF column). That panel is defined
in code but not yet built locally (only gnomAD v2.1.1 is materialized), so the first run builds
it; the 1000G EUR hg38 30x panel already on disk is an acceptable fallback if the v4.1 build is
impractically large. Source EAF: PPP
`A1FREQ` (frequency of ALLELE1 = EA); CSF `effect_allele_frequency`. One rule covers both
sources. Accepted cost: some palindromic SNVs and ambiguous indels near MAF 0.5 are dropped.

After both sides are in "NEA = hg38 REF at POS" form, trait and protein variants are matched
by an exact join on (CHR, POS, EA, NEA).

### Screening: a reusable HapMap3 protein-hit table, with a configurable threshold

Build once per proteomics source (not per trait) a table of every protein's significant
HapMap3 variants, from the existing slim per-protein HapMap3 files
(`HAPMAP_3_PPP_DATABASE`, `mecfs_bio/assets/gwas/ukbb_ppp/ppp_database/hapmap3/`) joined to
the shared HapMap3 index (`HAPMAP_3_PPP_DATABASE_INDEX`) for positions.

- **The p-value threshold is a parameter of this Task** (default 5e-8). Different thresholds
  are different assets.
- Store significant variants (protein id, CHR, POS hg38, p), not pre-merged loci, so a trait
  run can apply any window to it. (Builder may choose to also store a per-protein lead
  summary if convenient.)
- This replaces per-trait scans of ~2,940 files, which would be slow if the slim store lives
  on `/mnt/d`.

Note on cohort: the HapMap3 slim files on main are the **EUR discovery** cohort
(`eur_discovery_ppp_manifest.csv`, n ≈ 34,557). v1 uses discovery for both screening and the
dense download so they are consistent. (Memory notes say a Combined-cohort default was
introduced on a branch; if/when Combined slim files exist on main, the pipeline should take
the manifest/cohort as a parameter rather than hardcoding.)

### Colocalization method: coloc.abf first, behind an interface

v1 implements coloc.abf (Giambartolomei et al. 2014) in pure Python, vectorized over the
variants of a window. Defaults are coloc's: p1 = 1e-4, p2 = 1e-4, p12 = 1e-5. Priors are
parameters and are recorded in the raw table. Per-variant inputs are beta, se, and either
N + MAF (quantitative traits) or case fraction (case-control trait); the protein side is
always quantitative. The method interface takes the aligned window for one (trait, protein)
pair and returns posteriors; SuSiE-coloc later implements the same interface plus an LD
dependency.

The trait type (quantitative vs case-control) and, for case-control, the case fraction or
per-variant N_CASE/N_CONTROL must come from the trait meta / columns. coloc.abf's
case-control approximate Bayes factor uses the prior SD 0.2 on log-OR; quantitative uses
prior SD 0.15 × sdY. For proteins (inverse-rank normalized in PPP), sdY = 1.

## The pipeline (one trait run)

1. **Trait loci.** From the harmonized trait GWAS, find lead variants at p < 5e-8
   (parameter) and draw a window of ±500 kb (parameter) around each; merge overlapping
   windows. (Builder's choice whether to reuse `GwasLabLeadVariantsTask` or write a simple
   polars distance-based clumping; the latter avoids a gwaslab round-trip and is enough.)
2. **Candidate pairs.** Join trait windows with the protein-hit table: a protein is a
   candidate for a window if it has at least one significant variant inside it. Tag each pair
   **cis** if the protein's gene (hg38 coordinates from Sun et al. 2023 ST3, already parsed in
   `ppp_ldsc/gene_coords.py`) lies within ±1 Mb of the window, else **trans**. Trait windows
   overlapping the extended MHC (`extended_mhc_interval("38")`) are excluded by default
   (parameter) and recorded with status "excluded: MHC". The HapMap3 PPP index does NOT
   exclude the MHC on its own: it holds 2,149 variants in hg38 chr6:25.72–33.40 Mb (about half
   the density of the flanking 7.7 Mb), so MHC pQTL hits do reach the hit table.
3. **Select within a budget.** Score candidate pairs and pick which proteins to fetch (see
   "Candidate selection budget"). Every candidate not selected stays in the raw table with its
   score and status "not selected (budget)".
4. **Dense fetch, per selected protein.** Download the protein's raw tar from Synapse
   (`wf.download_from_synapse`), extract only the chromosome file(s) needed, slice every
   window that protein is a candidate for, then delete the tar and extracted files. A protein
   that is a candidate at several loci is downloaded once. (Optimization left to the builder:
   ranged-read only the needed chromosome member out of the tar via a presigned URL, as
   `fetch_synapse_file_head` already does for the first 64 KB; ~25 MB instead of ~545 MB. Not
   required for v1.)
5. **Orient.** Run each protein window through `resolve_chromosome` (untrusted path, gnomAD
   v4.1 panel); keep rows with null drop reason. Record flipped/dropped counts per pair.
6. **Align.** Exact join to the trait window on (CHR, POS, EA, NEA).
7. **Colocalize.** coloc.abf on the aligned window → H0–H4.
8. **Outputs.** Raw table, display table, plots (below).

## Candidate selection budget

Proposed by the user, refined together. High-N, polygenic quantitative traits can yield
hundreds of (locus, protein) candidates, most of them uninformative (hotspots, weak
signals). The run therefore tests only the highest-priority candidates, within a configured
budget.

- **The budget is a number of (trait locus, protein) pairs** (parameter, default 30). Each
  pair needs at most one fetch, so the pair count bounds worst-case cost for every source
  (PPP whole tar: 30 × ~545 MB ≈ 16 GB; PPP single-chromosome range read or CSF tabix: far
  less). A budget on distinct proteins was considered and rejected: it is only the right cost
  model when the fetcher downloads a whole protein file, and second loci are not free under
  range-read or tabix fetching, so it would leak one fetcher's cost model into the
  source-agnostic selection.
- **Fetch sharing is a fetcher optimization, not a selection rule.** When several selected
  pairs share a fetch unit (same protein tar; same protein + chromosome), the fetcher
  downloads it once.
- **Score per (trait locus, protein) pair**, from the trait GWAS and the hit table:
  - increases with the trait's −log10 p at the locus lead, **clipped** (default 50);
  - increases with the protein's min −log10 p inside the window, **clipped** (default 50),
    so p ≈ 1e-300 cis-pQTLs cannot swamp the score;
  - **pleiotropy penalty for trans pairs only**: decreases with the number of *other*
    proteins with hits in the window (from the hit table). Cis pairs are not penalized, so
    for example ABO's own cis-pQTL at the ABO locus stays a high-priority test. (Builder's
    choice: a cis bonus instead of exempting cis from the penalty, as long as cis pairs at
    hotspots are not de-prioritized.)
  - Exact functional form and weights are left to the builder, as parameters. A transparent
    additive form (e.g. clipped −log10 p's summed, minus w · log2(1 + n_other_proteins) for
    trans) is preferred over anything opaque.
- **Per-locus cap** (parameter, default 5 proteins per trait locus), so one strong locus
  cannot use up the whole budget.
- **Selection procedure**: rank pairs by score; walk down the ranking taking a pair while the
  pair budget and that locus's cap allow. Deterministic tie-breaking (e.g. by protein id), so
  the selection is reproducible.
- The score, its components and the selection status of every candidate go in the raw
  results table.

Memory: the window is the unit of work; at most one extracted chromosome file and one window
of panel rows are in memory at a time.

## Source interface (shape, not code)

Each proteomics source provides:

1. **Protein manifest**: stable protein id (PPP: OID, since control proteins share gene
   symbols across panels), gene symbol, and gene coordinates (for cis/trans).
2. **Screening store**: per-protein HapMap3 slim beta/se aligned to a shared variant index
   (exists for PPP; exists on the `csf-database` branch for CSF). The hit-table builder only
   reads this shape.
3. **Dense window fetcher**: (protein, list of windows) → dataframe per window in one common
   gwaslab-named schema — CHR, POS (GRCh38), EA, NEA, BETA, SE, EAF, N (per variant; constant
   for PPP) — declared genome build, untrusted. Orientation is NOT the fetcher's job; step 4
   does it uniformly.
4. **Sample-size semantics**: PPP N is constant per protein (present in the regenie `N`
   column); CSF N varies per variant. The common schema carries N per variant, so coloc.abf
   needs no special case.

## Outputs

### (a) Raw results table (parquet)

One row per (trait locus, protein) candidate pair: tested, not selected (budget), or failed
with a reason (e.g. too few aligned variants); plus one row per trait locus excluded as MHC.
Columns (builder may add):

- locus id, CHR, window start/end (hg38), trait lead variant (POS, EA, NEA, rsID if
  available) and its p;
- protein id, gene symbol, cis/trans, protein's lead variant in the window and its p;
- n variants in trait window, n in protein window, n aligned, n flipped / n dropped by
  orientation (by drop reason category);
- PP.H0–PP.H4, priors used (p1, p2, p12), method name;
- selection score and its components;
- status (tested / not selected (budget) / excluded: MHC / failed) and failure reason.

Integer columns Int32 (Tabulator cannot format parquet INT64).

### (b) Display table (parquet, docs-facing)

Rows from (a) with PP.H4 > 0.8 (parameter). Columns: protein (gene symbol, id), locus
(chr:start–end, lead variant), cis/trans, PP.H4. Written via
`dataframe_output.write_parquet_table` with zstd + byte-stream-split floats and Int32
integers, for `{{ data_table(...) }}`. An empty table (no hits) is a valid output, not an
error, so the docs page still builds.

### (c) Stacked locus plots

One figure per row of (b), PNG + SVG, matplotlib object-oriented API (no pyplot), modelled on
`polyfun_explain/polyfun_explain_plot_task.py`: panel 1 trait −log10 p, panel 2 protein
−log10 p, panel 3 genes via the existing `plot_gene_tracks` / `_plot_genes` path, shared
x-axis over the window. No LD colouring in v1. Highlight the trait and protein lead
variants. Title carries protein, cis/trans and PP.H4. Plots live in a result directory asset
alongside an index (protein id → filename) the docs can reference.

Which gene-info table to use (build-38 coordinates required) is the builder's choice among
existing gene-info assets; it must be GRCh38.

## Follow-ups (explicitly out of v1)

1. **CSF source.** Merge the needed parts of `csf-database` (manifest, HapMap3 index/slim
   builders). Dense fetch via the GWAS Catalog: either download-and-discard the ~200 MB
   `.tsv.gz`, or tabix-query the `harmonised/*.h.tsv.gz` (+ `.tbi`) for just the window over
   HTTP. Expect fewer hits (N ≈ 3,400).
2. **SuSiE-coloc.** One Broad UKBB LD window per trait locus (shared by every protein
   there), labels lifted 37→38; SuSiE on the trait once per locus and per protein; coloc on
   credible-set log-BFs. Enables LD r² colouring in plots.
3. **Combined PPP cohort** once its slim files exist on main.

## Scope, trait inputs, docs, resumability

- **v1 runs four traits**, to exercise the whole system and to make it likely that at
  least some runs produce real colocalizations (DecodeME may return none against PPP). All
  four have GRCh38 positions available (Willer LDL via the GWAS Catalog's harmonised
  liftover), so no liftover Task is needed in v1. (Liftover 37→38 is an acceptable
  preprocessing step for future traits.) For each, building its GRCh38
  `GenomeReferenceHarmonizationTask` asset (`UCSC_HG38_INDEXED_FASTA` +
  `GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES`) is part of this work.

  | Trait | Source asset | Type | Why included | Prep notes |
  |---|---|---|---|---|
  | ME/CFS (DecodeME) | `decode_me/raw_gwas_data`, build-38 form | case-control | the motivating trait | build-38 rsID asset exists; harmonize |
  | Seropositive RA (deCODE, Saevarsdottir 2022) | `DECODE_SEROPOSITIVE_RA_RAW` (`rheumtoid_arthritis/decode_seropositive/`) | case-control | immune trait; user's suggestion | `PosB38`; OR + P only, so BETA/SE via the existing `ComputeBetaPipe`/`ComputeSEPipe`; per-variant effective N via `EffectiveNFromCohortStringPipe` (as in `ra_seropositive_standard_analysis.py`); EAF is a percentage, use the `ScaleColPipe` step. Do NOT start from `SEROPOS_RA_FILTERED_FOR_FREQ` (MAF ≥ 0.05 drops variants coloc wants) |
  | LDL cholesterol (Willer et al. 2013, GLGC) | `mecfs_bio/assets/gwas/ldl/willer_et_al/raw/raw_willer_ldl_data.py` (GCST002222 harmonised `.h.tsv.gz`) | quantitative, per-variant N 50k–173k (median ≈ 90k) | polygenic trait with many loci: exercises the selection budget; lipid loci give likely positive controls. Chosen over MVP LDL, whose genetic correlations looked odd (possibly lipid-lowering drug use) | the GWAS Catalog file is GRCh38; use `hm_chrom`/`hm_pos` (7,641 of 2.43M rows failed liftover and are dropped). HapMap2-imputed, so sparser than the others. `hm_effect_allele_frequency` and `effect_allele_frequency` are 100% null; use `eaf_ref`, **verified** (`experiments/claude/coloc_allele_orientation/willer_eaf_ref_probe.py`, log `experiments/claude/logs/willer_eaf_ref_probe_hg38.log`) to be the 1000G EUR frequency of the ORIGINAL `effect_allele` (GLGC's Freq.A1.1000G.EUR): against the repo's 1000G EUR hg38 panel on 2.02M non-palindromic SNVs, corr 0.997, median |diff| 0.006, 99.8% within 0.05 (vs 4% for the other allele). It is a reference, not in-sample, frequency, which is adequate for coloc.abf and for panel-AF orientation. `beta` is relative to `effect_allele` (the Catalog's hm_beta sign flips exactly when its hm_effect_allele is the other allele). `base_pair_location` is also GRCh38 (equal to `hm_pos`, e.g. rs3094315 at chr1:817186), so either works. Recommended: take original `effect_allele`/`other_allele`/`beta`/`standard_error`/`eaf_ref`/`n` at `hm_chrom`/`hm_pos` and let our own harmonizer orient them (untrusted) rather than relying on the Catalog's `hm_*` alleles. Lowercase alleles must be uppercased |
  | Inflammatory bowel disease (Liu et al. 2023, EUR) | `LIU_ET_AL_2023_IBD_META_EUR_ONLY` | case-control | immune trait with many loci | raw is GRCh38 (the repo lifts it *to* 37 for other analyses; use the pre-liftover EUR-only asset); METAL-style alleles, so orientation is untrusted; AF_NFE, BETA_NFE, SE_NFE; no per-variant N, so use the existing prevalence/sample info in `liu_et_al_2023/auxiliary/` |

  Considered and rejected: CRP (Said et al. 2022) would be an obvious positive control (CRP
  is an Olink protein), but its file has no allele frequencies at all (both
  `effect_allele_frequency` and `hm_effect_allele_frequency` are 100% null over 10.9M rows),
  which both coloc.abf's quantitative path and untrusted harmonization need.

  **Expectation (Claude's guess, a sanity check rather than a requirement):** LDL should
  yield at least one high-H4 cis colocalization at a classic lipid locus with an Olink-measured
  protein (e.g. PCSK9). If LDL yields none, suspect the pipeline before the biology.
- **Trait type and case-control information** follow the existing pattern: a
  `BinaryPhenotypeSampleInfo`-style parameter (cf. `DECODE_ME_PREVALENCE_INFO` used by
  ct_ldsc) passed to the pipeline's `.create()`, absent for quantitative traits. No new field
  on `FilteredGWASDataMeta`. `.create()` asserts that the trait columns coloc.abf needs (BETA,
  SE, EAF, and N or the case-control info) are present.
- **Docs.** `docs/` is human-written and is not edited by this work. The deliverable stops at
  assets: the display-table parquet, the plot directory, and a plot index (protein id, locus
  id → filename). The handover includes suggested `data_table` / image macro calls for the
  user to paste into a page.
- **One Task per trait for the protein work.** The build graph is static (`deps` is fixed
  at graph construction), and which proteins to fetch depends on the contents of the hit
  table, so per-protein Tasks are not possible. All dense fetching, orienting and coloc for
  a trait run happen inside one Task. This is acceptable because the selection budget (see
  "Candidate selection budget") bounds the work: at most ~30 pairs × ~545 MB ≈ 16 GB, a few minutes
  of download, so a restart is cheap. Per-protein downloads use `call_with_retries`. (If a
  much larger budget is ever needed: a two-phase pattern, where the selection asset is
  exported to a committed CSV the way manifests are and per-protein Tasks are generated from
  it statically.)

## Left to the builder

- Task decomposition and asset layout (e.g. hit-table Task; per-trait candidate Task;
  per-trait coloc Task; display-table Task; plot Task), within repo conventions.
- Lead-variant/clumping implementation and window merge rule.
- Whether to implement the ranged-read single-chromosome fetch in v1.
- Concurrency of protein downloads (memory notes: ~4–8 concurrent is the useful range; the
  build-system scheduler is still serial, so in-Task concurrency or serial is fine).
- Handling of PPP "neither" SNVs that reverse-complement-match hg38 REF: `resolve_chromosome`
  strand-flips RC classes unconditionally (`resolve_chromosome.py:127`). PPP is plus-strand,
  so these are build drift, not strand errors; ≤ ~15k of 14.5M variants genome-wide. Either
  measure the effect or drop RC classes for sources known to be plus-strand.
- Minimum aligned-variant count for a pair to be tested.
- Choice of gnomAD ancestry column (NFE expected) via the panel's `panel_ancestry`.

## Calls made by Claude (not explicitly decided by the user)

- ±500 kb trait windows; ±1 Mb gene distance for cis.
- p < 5e-8 for trait leads (parameter).
- PP.H4 > 0.8 for the display table (parameter).
- coloc default priors (p1 = p2 = 1e-4, p12 = 1e-5), as parameters.
- Hit table stores significant variants rather than merged loci.
- EUR discovery cohort for v1.
- gnomAD v4.1 genomes (GRCh38) as the orientation panel.
- v1 traits besides DecodeME: seropositive RA and Willer LDL (both user's choice), Liu 2023
  IBD (Claude's); CRP rejected for missing allele frequencies.
- Case-control info passed as a BinaryPhenotypeSampleInfo-style parameter, not a new meta
  field.
- Deliverable stops at assets; no edits to docs/.
- Selection defaults: 30 pairs, 5 per locus, −log10 p clipped at 50; MHC excluded by
  default. (The user approved the scheme and the refinements; the default numbers are
  Claude's.)
