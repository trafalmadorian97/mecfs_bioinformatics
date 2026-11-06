# CSF protein colocalization (Western et al. 2024) — design

**Date:** 2026-11-05
**Status:** approved 2026-11-06 (builder check done 2026-11-06)

## Intent

Main now has large-scale coloc.abf between a trait GWAS and every UKB-PPP plasma protein with a
pQTL near a trait locus (`mecfs_bio/build_system/task/protein_coloc`, merged in #1327). We want the
same analysis against the Western et al. 2024 CSF proteomics GWAS (SomaScan 7k, 7,008 aptamers,
n ≈ 3,506 European ancestry, GRCh38; Nature Genetics 56:2672, PMID 39528825).

Why: the unmerged `csf-database` branch built a CSF pQTL database for trait × CSF-protein genetic
correlation, but n ≈ 3.5k gave no significant genetic correlations, so it was never merged. Coloc
depends on strong cis-pQTL signal rather than genome-wide polygenic signal, and CSF cis-pQTLs
should still be strong at that sample size. CSF is the compartment closest to brain, so
neurological / psychiatric traits are where this should add something PPP (plasma) cannot.

Success: for each trait, the same three outputs the PPP run produces (raw results, PP.H4 display
table, stacked locus plots), runnable with one call per trait from a driver module.

**Alzheimer's is a positive-control gate (user decision).** Western et al. framed their paper
around AD causal proteins. AD runs first, and its results are checked against the AD
proteins Western et al. report (their coloc / MR results — pull the list from the paper and its
supplementary tables; read the PDF, per repo practice, don't work from memory or the abstract).
Commit the extracted list as a small reference file with its provenance. The check compares, per
reported protein: was it a candidate, was it selected within the pair budget, did it reach
PP.H4 ≥ 0.8. Misses not explained by method differences (Western's own coloc settings, windows,
priors) or by the budget are investigated as possible bugs (fetcher, protein ids, orientation)
before the other seven traits run. The outcome of the check is reported to the user; the
builder sets the exact comparison script (in experiments/claude) but not a pass threshold — the
user judges it.

## Goals

1. Bring the **database part** of the `csf-database` branch onto main, unchanged in asset ids and
   on-disk paths, so the 7,008 already-built slim HapMap3 aptamer files are reused, not rebuilt.
2. **Generalize** the protein-coloc machinery so PPP and CSF share classes wherever the logic is
   source-independent (hit table, coloc task, asset-generator group, protein id typing).
3. Add a **CSF dense-window fetcher** that reads trait-locus windows remotely via tabix.
4. Add a **CSF hit table** asset and a **CSF coloc asset generator**, mirroring the PPP ones.
5. Run against eight traits:
   - Existing PPP traits (have GRCh38-harmonized tasks): DecodeME GWAS-1, deCODE seropositive RA,
     Willer LDL, Liu et al. 2023 IBD (EUR).
   - Existing in repo but need a GRCh38 genome-reference-harmonized task: Alzheimer's
     (Bellenguez et al. 2022), schizophrenia (PGC 2022, Trubetskoy et al.).
   - New ingestion: Parkinson's (Nalls et al. 2019, without 23andMe) and ALS (van Rheenen et al.
     2021, EUR only).
6. A `mecfs_bio/analysis/csf_coloc_analysis.py` driver with one `run_*` function per trait, like
   `ppp_coloc_analysis.py`.
7. **Exploratory DecodeME run with a relaxed locus threshold** (user request): in addition to the
   standard DecodeME × CSF run at lead p < 5e-8, a second DecodeME × CSF run whose
   `ColocRunConfig.lead_p_threshold` is relaxed so sub-genome-wide-significant ("tentative")
   DecodeME loci are also tested against CSF pQTLs. Separate asset ids (e.g. base name
   `decode_me_gwas_1_lead_p_<threshold>`), same three outputs, its own `run_*` in the driver, and
   clearly labelled exploratory in the module docstring. The threshold value is the user's to set
   (user decided: **two exploratory runs, lead p < 1e-6 and lead p < 1e-5**, each with its own
   asset ids and driver function); the pair budget in `SelectionConfig` may need raising to match the larger
   number of loci.

## Non-goals

- The CSF LDSC heritability / genetic-correlation code on `csf-database` (`csf_ldsc/`,
  `batched_ldsc/` refactor, `hapmap3_csf_heritability.py`, its experiments). Stays on the branch.
- Any text under `docs/`. Docs pages are human-written; the user will reference the outputs.
- Correcting for sample overlap between Western's cohorts and the trait GWASs (see Risks).
- SuSiE-coloc or any multi-signal method; this is coloc.abf only, as for PPP.

## Anti-goals (would make this a failure even if it runs)

- Rebuilding the PPP coloc outputs or the CSF slim database because asset ids / paths changed in
  the generalization. PPP asset ids, metas and on-disk paths must be byte-identical after the
  refactor; CSF slim asset ids and paths must match what the branch produced.
- Bulk-downloading the ~1.4 TB of dense CSF summary statistics.
- Changing PPP coloc numerical results.

## Constraints

- Coding conventions in the user's memory/CLAUDE.md apply (attrs frozen, Literal/NewType ids,
  no monkeypatch — inject fetchers, imports at top, Task-level tests, docstrings without
  backticks, `pixi r invoke green` after significant changes).
- 16 GB RAM machine; dense data must be windowed.
- The verifying-trace build system hashes (imohash: three 16 KB samples per file) every file on
  a requested output's dependency graph. The CSF hit table depends on all 7,008 slim files, so
  they must live on the local disk, not the /mnt/d external drive (per-file latency there makes
  the trace check unacceptably slow — user's decision).

## Decided by the user

- Alzheimer's is a gate checked against Western et al.'s published AD proteins.
- Exploratory DecodeME × CSF runs at lead p < 1e-6 and < 1e-5, in addition to 5e-8.
- Trait set above; Parkinson's source is **Nalls 2019 (GCST009325)**, not Kim 2024.
- ALS source: van Rheenen 2021 **EUR-only (GCST90027164)**.
- New GWAS ingestion tasks are in scope for this project.
- Bring over only the database part of `csf-database`.
- Move the CSF slim files back from /mnt/d to the main disk.
- Same outputs as PPP: raw results + display table + locus plots per trait.
- Generalize shared pieces so PPP and CSF share classes where sensible.

## Facts established during brainstorming (verified 2026-11-05)

### CSF database already built

- `default_runner_config.yaml` (gitignored, machine-local) has a `path_remap` rule sending
  `gwas/western_csf_pqtl` to `/mnt/d/asset_store_remote`. All **7,008** slim files exist there
  (`gwas/western_csf_pqtl/<GENE>/aligned/hapmap_3_index/csf_slim_hapmap_3_<GENE>_<seqid>.parquet`,
  ~6.5 MB each, ~45 GB, 6,179 gene folders). An old `csf_heritability/` output is also there.
- The CSF variant index (`reference_data/csf_pqtl_variant_index/hapmap_3_membership_list/processed/csf_variant_index.parquet`,
  12.5 MB) and the template aptamer (`reference_data/csf_pqtl/western_2024/...`) are on the
  local disk already.
- Branch files constituting "the database part":
  - `mecfs_bio/build_system/task/csf_database/{build_slim_aptamer_parquet_task,construct_csf_variant_index_task,gwas_catalog_url}.py`
  - `mecfs_bio/asset_generator/csf_slim_aptamer_asset_generator.py`
  - `mecfs_bio/assets/gwas/csf_pqtl/raw/csf_template_aptamer.py`
  - `mecfs_bio/assets/gwas/csf_pqtl/csf_database/hapmap3/{hapmap3_csf_index,hapmap3_csf_database_aptamer_files}.py`
  - `mecfs_bio/assets/reference_data/csf_pqtl_sumstats/{csf_aptamer_manifest.csv,regenerate_csf_manifest.py}`
  - `mecfs_bio/constants/csf_database_constants.py`
  - `mecfs_bio/analysis/csf_pqtl_database_build.py`
  - the matching tests under `test_mecfs_bio/unit/build_system/task/csf_database/`
  - any helper these import that is not on main (the branch also touched
    `ppp_database/allele_key.py`); bring over only what the database files need. Known so far
    (from the builder check): `gwas_ssf_constants.py`, and the branch code must be adapted to
    main's current conventions (e.g. `unordered_allele_key` now expects an `assert_all_snv`
    guard; the `request_connections` parameter added to WF and the downloaders on main).
  - If the database does not verify cleanly (e.g. because an upstream asset such as
    `HAPMAP_3_MEMBERSHIP_LIST` was rebuilt since the branch ran), **stop and report** to the
    user before letting anything rebuild — a full CSF rebuild is ~1.4 TB of downloads (my call).
    The same rule applies if PPP coloc outputs fail to verify after the refactor.
  - The plan that built it: `experiments/claude/design_specs/2026-08-19-csf-pqtl-database-plan.md`
    on the branch.
- Manifest columns: analyte, seq_id, uniprot, entrez_gene_symbol, target_full_name, accession, md5.
  78 aptamers have a multi-gene symbol.

### Moving the files back

Automated: `pixi r invoke migrate-asset-store --reverse --only western_csf_pqtl` (dry run), then
with `--execute`. It rsyncs, verifies, deletes the /mnt/d copy. Only **after** that, remove the
`gwas/western_csf_pqtl` line from `default_runner_config.yaml` (the task only knows prefixes the
config still declares). This is an operational step the user runs; it is not code.

### CSF dense data via remote tabix

- Every Western study on the GWAS Catalog has `harmonised/<GCST>.h.tsv.gz` plus a `.tbi`; the
  FTP server supports byte ranges. URL:
  `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/<bucket>/<GCST>/harmonised/<GCST>.h.tsv.gz`
  (`gwas_catalog_bucket` from the branch derives the bucket).
- `pysam.TabixFile(url).fetch("19", 44_400_000, 45_400_000)` returned 2,808 rows in ~4 s
  (pysam is in the default env).
- Build: GRCh38, 1-based (verified from the study meta yaml: `genome_assembly: GRCh38`).
- Header: `chromosome base_pair_location effect_allele other_allele beta standard_error
  effect_allele_frequency p_value variant_id rsid neg_log_10_p_value rs_id n
  hm_coordinate_conversion hm_code`. Contigs are "1".."22" (no chrX, no "chr" prefix).
  GRCh38. Per-variant `n` (~3,150–3,450).
- ~1.6% of variants failed GWAS Catalog harmonisation and are absent from the harmonised file
  (present in the raw, unindexed file). Accepted (my call).
- The harmonised effect allele may differ from the raw file's; irrelevant because protein windows
  go through the existing untrusted genome-reference orientation anyway. (Note from earlier work:
  the raw file's effect allele is the in-sample minor allele.)

### New trait sources (GWAS Catalog, verified)

| Trait | Accession | Sample | Build | Key columns | Size |
|---|---|---|---|---|---|
| Parkinson's (Nalls 2019, ex-23andMe) | GCST009325 | 15,056 cases + 18,618 UK Biobank proxy cases, 449,056 controls, EUR | GRCh37 | chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, **N_cases, N_controls** (per variant) | 0.96 GB `.tsv` (uncompressed) |
| ALS (van Rheenen 2021, EUR) | GCST90027164 | 27,205 cases, 110,881 controls | GRCh37 | rsid, chromosome, base_pair_location, effect_allele, other_allele (**lowercase**), effect_allele_frequency, beta, standard_error, p_value, **N_effective** (per variant) | 197 MB `_buildGRCh37.tsv.gz` |

Raw file URLs:
- `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/GCST009325.tsv` (md5 in `md5sum.txt` beside it; meta yaml `data_file_md5sum: 88c9f0d63ebeda700396d4402d8adf03`)
- `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz` (md5 in `md5sum.txt`)

Existing traits needing a GRCh38 harmonized task:
- Alzheimer's: `BELLENGUEZ_ET_AL_ALZHEIMERS_RAW` downloads the GWAS Catalog *harmonised*
  `GCST90027158` file (already GRCh38 hm_ coordinates); `BELLINGUEZ_ET_AL_ALZHEIMERS_EXTRACTED`
  unzips it.
- Schizophrenia: `PGC_2022_SCH_RAW` (figshare, `PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz`,
  GRCh37).
- Pattern to follow: `ra_seropositive_build_38_harmonized.py`, `willer_ldl_build_38_harmonized.py`,
  `liu_et_al_2023_eur_build_38_harmonized.py`, DecodeME's `decode_me_gwas_1_build_38_harmonized.py`
  (gwaslab sumstats → table → `GenomeReferenceHarmonizationTask` against UCSC hg38 FASTA and the
  gnomAD v4.1 nfe panel). Keep indels; filter MAF ≥ 0.01 only if needed for memory, as RA did.
- Correction from the builder check: there is **no existing 37 → 38 liftover** pattern. RA, IBD
  and DecodeME are GRCh38 at source; Willer LDL uses the GWAS Catalog harmonised file's GRCh38
  `hm_chrom`/`hm_pos`. The only repo liftovers go 38 → 37. Route for GRCh37 sources (my call):
  use the GWAS Catalog harmonised `.h.tsv.gz` (GRCh38) as Willer does where one exists — it does
  for Nalls (`GCST009325/harmonised/GCST009325.h.tsv.gz`) and van Rheenen
  (`GCST90027164/harmonised/34873335-GCST90027164-MONDO_0004976.h.tsv.gz`) and is already what
  the Bellenguez download is. The builder must confirm the harmonised files keep the per-variant
  N columns (N_cases/N_controls, N_effective); if not, join them back from the raw file. For
  schizophrenia (figshare, GRCh37) the user believes no native GRCh38 version exists: add a
  gwaslab 37 → 38 liftover step (user approved; also fine for any other trait that needs it).

## The how

### 1. Database (from the branch)

Cherry-pick / copy the files listed above. Asset ids, metas, groups and paths must be exactly as
on the branch so the existing files and their traces in the info store are recognized. Verify by
running the build for the CSF slim collection after the move and confirming nothing
re-materializes (dry/verify-only run, or a run on a handful of aptamers with logging).

### 2. Generalization of protein_coloc (shared by PPP and CSF)

Source-independent pieces become shared classes with neutral names. Indicative (builder decides
exact names):

- **Protein id**: today `Oid` (Olink id) flows through candidate selection and gene coords.
  Introduce a neutral protein-id type (e.g. a NewType over str) for the coloc layer; PPP passes
  Oids, CSF passes analytes/seq ids.
- **Hit table task**: `PppProteinHitTableTask` takes PPP slim tasks. Generalize to take any
  collection of slim per-protein tasks over a shared variant index plus a protein id / gene
  label per task. PPP and CSF both instantiate it.
- **Gene coordinates**: the coloc task currently reads Sun et al. ST3 via `read_gene_coords`.
  Generalize to a mapping protein id → GRCh38 gene interval, supplied by a source-specific reader:
  PPP keeps ST3; CSF maps the manifest's `entrez_gene_symbol` through the GRCh38 MAGMA gene
  location file (`MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED`). Multi-gene or
  unmatched aptamers get no coordinates → cis status "unknown", exactly as PPP complexes do (my
  call).
- **Coloc task**: `TraitProteinColocTask` already depends on the `DenseWindowFetcher` protocol;
  keep it, plus the generalized hit table / gene coords, and make the protein trait model a
  parameter (both sources use sd_y = 1).
- **Asset generator**: `generate_ppp_coloc_tasks` / `PppColocTaskGroup` become a shared generator
  + group parameterized by hit table, fetcher, gene-coords source and an asset-id infix
  (`ppp` / `csf`); thin `generate_ppp_coloc_tasks` / `generate_csf_coloc_tasks` wrappers default
  the source-specific inputs.
- **Display table and locus plots**: should work unchanged once ids are neutral; the gene label
  for CSF is the manifest gene symbol (aptamer id alongside, since one gene can have several
  aptamers).
- **Hard rule**: PPP asset ids (`{base}_ppp_coloc`, `_display`, `_plots`, the hit table's id),
  metas and paths unchanged. Existing PPP tests keep passing.

### 3. CSF dense-window fetcher

`CsfWindowFetcher` implementing `DenseWindowFetcher`: for one aptamer, open the remote harmonised
file with pysam TabixFile (URL from the manifest accession), fetch each requested window, map
columns to `WINDOW_SCHEMA` (CHR, POS, EA, NEA, BETA, SE, EAF, N), return frames. Retries via
`call_with_retries` as PPP does. Chromosome 23 windows return empty (no chrX in CSF). No local
download. Unit tests use a small local bgzipped+tabixed fixture (pysam accepts a path), injected;
no network in unit tests.

### 4. Protein model

Western z-score-normalized aptamer levels → `QuantitativeTraitModel(sd_y=1.0)`, same as PPP.
N per variant from the `n` column.

### 5. CSF hit table and generator wiring

- `mecfs_bio/assets/gwas/csf_pqtl/csf_database/hapmap3/hapmap3_csf_hit_table.py`: the
  generalized hit table over the 7,008 CSF slim tasks, p < 5e-8 (same threshold as PPP; my call).
- `mecfs_bio/asset_generator/csf_coloc_asset_generator.py` (or the shared generator's CSF
  wrapper) defaulting hit table, `CsfWindowFetcher` from the manifest, CSF gene coords, UCSC hg38
  FASTA, gnomAD v4.1 nfe panel, MAGMA gene track.
- Per-trait modules `.../analysis/<trait>_csf_coloc.py` next to each `<trait>_ppp_coloc.py`
  (new ones for AD, SCZ, PD, ALS).
- `ColocRunConfig` defaults (±500 kb windows, MHC excluded, min 50 aligned variants) as PPP.
  The PPP pair budget (`SelectionConfig.max_pairs = 30`, `max_per_locus = 5`) exists because
  each PPP pair costs a full Synapse tar download; a CSF pair costs one remote tabix read
  (~4 s), so the CSF generator defaults to a larger budget (user agreed; the value is the
  builder's call, stated with its reasoning in the generator docstring), and the exploratory
  DecodeME runs (goal 7) need one sized to their loci.

### 6. Trait models

- DecodeME, RA, LDL, IBD: same trait model as their PPP modules.
- Alzheimer's, schizophrenia: case-control with the GWAS's case fraction (builder derives from
  the sources; per-variant N where provided).
- Parkinson's: count proxy cases as cases, case_fraction = (15,056 + 18,618) / 482,730 (my call;
  proxies dilute effects slightly, effect on PP4 expected small). Per-variant N = N_cases +
  N_controls.
- ALS: `N_effective` as N with case_fraction = 0.5 (the convention effective N assumes) (my call).

## Left to the builder

- Exact class / module names of the generalized pieces and how far the rename reaches, within
  the "PPP ids/paths unchanged" rule.
- Whether the Bellenguez harmonised GRCh38 coordinates are used directly or re-derived.
- How SCZ / PD / ALS get to GRCh38 (follow the existing b37 → b38 harmonized trait pattern).
- Memory-saving MAF filters for large trait files, if needed on 16 GB.
- Test design (Task-level, injected fetchers, small fixtures).

## Risks / caveats to record in outputs or docstrings

- **Sample overlap**: Western's CSF cohorts (Knight ADRC, ADNI, …) likely overlap Bellenguez and
  possibly Nalls. coloc.abf assumes independent samples; not corrected.
- **Power**: n ≈ 3.5k means mostly cis-pQTLs will pass the hit-table screen; few trans pairs.
- **Multiple aptamers per gene**: results are per aptamer; display shows gene + aptamer id.
- **Missing 1.6%**: harmonised-file variants that failed Catalog harmonisation are absent from
  protein windows.

## Operational order

1. User runs `migrate-asset-store --reverse --only western_csf_pqtl --execute`, then deletes the
   path_remap line.
2. Database code onto main; confirm the slim collection verifies without rebuilding.
3. Generalization refactor; confirm PPP tests unchanged and PPP outputs verify without rebuild.
4. CSF fetcher, gene coords, hit table, generator.
5. Trait ingestion / harmonization for AD, SCZ, PD, ALS.
6. Run Alzheimer's; positive-control check against Western et al.'s AD proteins; report to the
   user.
7. Run the other seven traits and the two exploratory DecodeME runs.
