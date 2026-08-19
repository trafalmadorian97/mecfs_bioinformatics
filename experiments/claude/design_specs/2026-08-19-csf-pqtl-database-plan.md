# CSF pQTL Database — Implementation Plan

**Date:** 2026-08-19
**Status:** Awaiting review, not yet implemented
**Scope:** Build the per-aptamer summary-statistics database for Western et al. 2024 CSF
proteogenomics, analogous to the existing UKB-PPP database under
`mecfs_bio/build_system/task/ppp_database/`.

**Out of scope (deliberately deferred):** LDSC heritability / rg over the CSF database. That
work depends on an unresolved statistical question (see [Deferred decisions](#deferred-decisions))
and is not needed to build the database.

---

## Goal

Produce, for each of 7,008 CSF aptamers, a slim parquet file holding `BETA`, `SE` and `N` as
`float32` in the row order of a shared variant index — the same storage contract as the PPP
database, so the same downstream machinery can consume it later.

---

## Source facts (verified 2026-07-27 / 2026-08-19)

Everything below was measured, not inferred. It is restated here so this plan stands alone.

### Where the data lives

Two sources exist. **Use the GWAS Catalog.**

| | Box (linked in the paper) | GWAS Catalog |
|---|---|---|
| Contents | 7,009 raw PLINK2 `X<SeqId>.glm.linear.gz`, ~225 MB each, 1.58 TB | 7,008 studies, GWAS-SSF v1.0, ~203 MB each |
| Listing | scrape the `Box.postStreamData` JSON blob, 351 pages × 20 | REST API |
| Checksums | none published | `md5sum.txt` + `data_file_md5sum` in meta.yaml |
| Transport | anonymous HTTP | plain HTTPS — `wf.download_from_url` works unchanged |

Paper: Western et al., *Nature Genetics* 56:2672 (2024), PMID **39528825**. 3,506 unrelated
European-ancestry individuals across 8 cohorts, SomaScan 7k, GRCh38, TOPMed r2 imputation,
PLINK2 linear regression on z-score-normalised aptamer levels with age, sex, 10 PCs, cohort
and array as covariates.

URL pattern (accessions run `GCST90421001`–`GCST90428xxx`, so the bucket segment varies):

```
https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/
    GCST{lower}-GCST{upper}/{accession}/{accession}.tsv.gz
```

The bucket is the 1000-wide range containing the accession, e.g. `GCST90421540` lives under
`GCST90421001-GCST90422000/`. Derive it arithmetically; do not hardcode a list.

### File shape

GWAS-SSF v1.0, tab-separated, gzipped. Columns:

```
chromosome  base_pair_location  effect_allele  other_allele  beta
standard_error  effect_allele_frequency  neg_log_10_p_value  variant_id  rs_id  n
```

Measured on `GCST90421540` (7,327,953 variants):

| Property | Value |
|---|---|
| Chromosomes | **autosomes 1–22 only; no chrX** |
| Build | GRCh38, 1-based |
| Sort order | chromosome, then position — already sorted |
| Multiallelic sites | **none** — `(chrom, pos)` is unique |
| `rs_id` | complete |
| Indels | 368,652 |
| MAF floor | ~0.001 |
| `n` | **varies per variant**: 3,153–3,446, median 3,446, 230 distinct values |
| HapMap3 overlap | **1,025,155 / 1,206,640 = 85.0%** (allele-aware) |

**Neither the observed sort order nor the observed `(chrom, pos)` uniqueness may be relied
on.** Both are properties of one file that we happen to have inspected, not guarantees of the
format, and the index row order is the alignment slot for every aptamer in the database — a
silent reordering corrupts all 7,008 files at once. Always sort explicitly on the full key
`(chromosome, position, effect_allele, other_allele)`, and never on `(chromosome, position)`
alone, even though that pair is unique today.

Note the HapMap3 coverage is far below PPP's 99.24%, so the CSF index is a genuinely
different variant set. Any future joint PPP × CSF analysis needs an intersection index; do
not assume the two indices are interchangeable.

Because GWAS-SSF names the effect allele explicitly (`effect_allele` / `other_allele`), the
PLINK2 `A1` vs `REF`/`ALT` orientation question that applies to the Box files does not arise
on this route.

### Aptamer identity, and the accession → aptamer map

`aptamer_info.xlsx` (826 KB, at the Box link root) has 7,584 rows × 19 columns; the relevant
ones are `Analytes` (e.g. `X10000.28`), `SeqId` (`10000-28`), `UniProt`, `EntrezGeneSymbol`,
`TargetFullName`, `Organism`, `Type`, `Step Removed`. The 576 rows with a non-null
`Step Removed` are exactly the ones not published: **7,584 − 576 = 7,008**.

`Analytes` is the primary key. Gene symbol is **not** unique — 7,584 aptamers cover only
6,398 gene symbols — so asset ids must carry the SeqId, mirroring PPP's `gene_OID` pattern.

The Catalog assigns one accession per aptamer with a **unique** trait string, having already
disambiguated shared target names itself. A four-rule resolver yields a complete bijection:

| Rule | Resolves |
|---|---|
| 1. Parse `(analyte X####.##)` out of the trait string | 1,256 |
| 2. Exact `"<TargetFullName> levels"` | 5,739 |
| 3. Case-insensitive fallback (`Thyroxine-Binding Globulin` vs `-binding`) | 10 |
| 4. Hardcoded override, 3 Casein kinase II naming variants | 3 |
| | **7,008 / 7,008** |

The three overrides, which must be a module-level constant in the manifest generator:

| Catalog trait | `Analytes` | Gene(s) |
|---|---|---|
| `Casein kinase II subunit alpha-2 levels` | `X13681.173` | CSNK2A2 |
| `Casein kinase II alpha-1: beta heterotetramer levels` | `X5225.50` | CSNK2A1\|CSNK2B |
| `Casein kinase II alpha-2: beta heterotetramer levels` | `X5226.36` | CSNK2A2\|CSNK2B |

### Storage, measured

Probe: `experiments/claude/csf_pqtl_database/n_storage_cost_probe.py`, log at
`experiments/claude/logs/csf_n_storage_cost_probe.log`. Aligned to HapMap3 (1,025,155 rows).

Cost of storing per-variant `N`, at zstd level 22:

| Encoding | File | `n` column | vs beta+se | × 7,008 |
|---|---|---|---|---|
| beta, se | 6,646,592 B | — | — | 46.6 GB |
| **beta, se, N — all byte-stream-split** | **6,766,058 B** | **119,319 B** | **+1.8%** | **47.4 GB** |
| beta, se BSS + N float32 dictionary | 6,790,551 B | 143,800 B | +2.2% | 47.6 GB |
| beta, se BSS + N int32 dictionary | 6,790,521 B | 143,750 B | +2.2% | 47.6 GB |

Per-variant `N` costs **+1.8%**, i.e. 0.8 GB across the whole database. Store it
unconditionally. Byte-stream-split beats dictionary encoding on the aligned subset (the
reverse of the full 7.3M-row file, where dictionary wins 1.74 MB vs 2.19 MB) — once sorted
and thinned, `N`'s high bytes are near-constant and the byte planes compress extremely well.
**Do not re-derive this from full-file numbers; they point the wrong way.**

Compression level, same shape:

| Level | Size | Write | × 7,008 |
|---|---|---|---|
| default | 7.22 MB | 0.10 s | 50.6 GB, 0.2 h |
| 9 | 7.08 MB | 0.17 s | 49.6 GB, 0.3 h |
| **15** | **6.78 MB** | **1.08 s** | **47.5 GB, 2.1 h** |
| 19 | 6.77 MB | 2.69 s | 47.4 GB, 5.2 h |
| 22 | 6.77 MB | 4.22 s | 47.4 GB, 8.2 h |

**Use level 15.** It captures essentially all the available gain; 19 and 22 spend 3–6 extra
CPU-hours across the database to save ~16 MB total.

Level does **not** affect downstream read speed, which is the property that matters for LDSC.
Reading all columns into numpy (LDSC's access pattern), 25 reps, median with p10–p90:

| Level | Size | Cold read | Warm read |
|---|---|---|---|
| default | 7.22 MB | 23.9 ms (18.7–33.1) | 17.2 ms (15.2–29.1) |
| 3 | 7.17 MB | 17.4 ms (16.2–23.7) | 14.0 ms (12.4–18.8) |
| 9 | 7.08 MB | 15.7 ms (15.2–17.4) | 12.1 ms (11.6–14.3) |
| **15** | **6.78 MB** | **16.9 ms (16.2–20.6)** | **13.2 ms (12.7–15.8)** |
| 19 | 6.77 MB | 17.0 ms (16.0–27.9) | 13.7 ms (12.7–17.9) |
| 22 | 6.77 MB | 17.5 ms (16.3–22.4) | 14.7 ms (13.0–20.5) |

The spread across levels 3–22 (~1–2 ms) is no larger than the run-to-run spread within a
single level, so the levels are indistinguishable. This is the expected zstd property:
decompression speed is largely independent of compression level, and a higher level gives a
smaller file to read. At ~17 ms per aptamer, one pass over all 7,008 slim files costs about
two minutes — negligible against the regression itself.

---

## Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Source | GWAS Catalog | Stable URLs, md5 checksums, standard format, no HTML scraping, no new `WF` capability |
| Existing PPP tasks | **Not modified** | User constraint. Also: code changes do not invalidate the build cache, so refactoring the PPP index writer would silently leave a behaviour change sitting in the 2,944 built slim files |
| Code sharing | New sibling package, some duplication | ~100 lines overlap with `ConstructPppVariantIndexTask`. Unifying later is cheap and can be done deliberately with a forced rebuild |
| Primary key | `Analytes` (`X####.##`) | Gene symbol is not unique across aptamers |
| Manifest | Committed CSV + offline regenerator | Matches `regenerate_ppp_manifest.py`. Catalog REST paging takes ~13 min, too slow for build time |
| `N` | Stored per variant, `float32`, byte-stream-split | +1.8% storage; keeps the statistical decision open at zero cost |
| Compression | zstd level 15, BSS on all three float columns | Measured sweet spot |
| Index membership | Reuse `Hapmap3MembershipTask` unchanged | Already produces the `(CHR, POS, EA, NEA, rsID)` contract |
| hg19 position | **Omitted from the index** | GWAS-SSF carries only GRCh38; PPP got hg19 free from the regenie `ID`. Add a liftover task later only if something needs it |

---

## Architecture

### Task graph

```
DownloadFileTask (template aptamer, GCST90421540.tsv.gz)
        |
        v
CompressedCSVToParquetTask  ------------+
                                        |
Hapmap3MembershipTask (existing) -------+--> ConstructCsfVariantIndexTask
                                                       |
                                                       v
                                        BuildSlimCsfAptamerParquetTask  (x 7,008)
                                        (downloads its own .tsv.gz, no asset dep on it)
```

Each slim task downloads its aptamer into its scratch dir and discards it, exactly as
`BuildSlimProteinParquetTask` does — the 203 MB source is never materialised as an asset.

### New files

| Path | Purpose |
|---|---|
| `mecfs_bio/constants/gwas_ssf_constants.py` | GWAS-SSF v1.0 column names (a public standard, not CSF-specific) |
| `mecfs_bio/constants/csf_database_constants.py` | `Analyte`/`SeqId`/`GcstAccession` NewTypes, CSF index column names |
| `mecfs_bio/build_system/task/csf_database/__init__.py` | |
| `mecfs_bio/build_system/task/csf_database/construct_csf_variant_index_task.py` | `ConstructCsfVariantIndexTask` |
| `mecfs_bio/build_system/task/csf_database/build_slim_aptamer_parquet_task.py` | `BuildSlimCsfAptamerParquetTask`, `CsfAptamerFile` |
| `mecfs_bio/build_system/task/csf_database/gwas_catalog_url.py` | Accession → FTP URL (bucket arithmetic) |
| `mecfs_bio/asset_generator/csf_slim_aptamer_asset_generator.py` | One task per manifest row |
| `mecfs_bio/assets/reference_data/csf_pqtl_sumstats/regenerate_csf_manifest.py` | Offline manifest builder |
| `mecfs_bio/assets/reference_data/csf_pqtl_sumstats/csf_aptamer_manifest.csv` | Committed, 7,008 rows |
| `mecfs_bio/assets/gwas/csf_pqtl/raw/csf_template_aptamer.py` | Download + parquet-convert the template |
| `mecfs_bio/assets/gwas/csf_pqtl/csf_database/hapmap3/hapmap3_csf_index.py` | The index asset |
| `mecfs_bio/assets/gwas/csf_pqtl/csf_database/hapmap3/hapmap3_csf_database_aptamer_files.py` | Terminal collection |
| `mecfs_bio/analysis/csf_pqtl_database_build.py` | `DEFAULT_RUNNER.run(...)` entry point |
| `test_mecfs_bio/unit/build_system/task/csf_database/test_*.py` | Unit tests |

Manifest columns: `analyte, seq_id, uniprot, entrez_gene_symbol, target_full_name, accession`.

---

## Implementation steps

Each step ends green (`pixi r invoke green`, output tee'd to a logfile) before the next
starts.

### Step 1 — Constants

`gwas_ssf_constants.py` and `csf_database_constants.py`. No behaviour; unblocks everything
else and keeps column-name literals out of the task bodies.

### Step 2 — Manifest

`regenerate_csf_manifest.py`: page the Catalog REST API for PMID 39528825, download
`aptamer_info.xlsx` from the Box link, apply the four resolver rules, write the CSV.

Assertions (all fatal, per the fail-fast convention):
- exactly 7,008 accessions returned;
- exactly 7,008 aptamers survive the `Step Removed` filter;
- the resolver terminates at a **complete bijection** — every accession maps, every kept
  aptamer is hit, no analyte used twice.

Run once (~15 min); commit the CSV. Re-run only to pick up upstream re-curation.

### Step 3 — Template aptamer asset

`DownloadFileTask` for `GCST90421540.tsv.gz`
(md5 `047befd46b553da2bcecf7c8faa91749`), then `PipeDataFrameTask` to convert it to parquet.
Mirrors `STACK_UKBBPPP_RABGAP1L`'s role for PPP.

`PipeDataFrameTask` reads the whole 7.3M-row file into memory, which is acceptable for this
one-off template asset. Two constraints follow from its implementation:

- the `DownloadFileTask` meta must carry
  `read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t"))`, since
  `PipeDataFrameTask` reads its source through `scan_dataframe_asset`;
- `backend="polars"` is required — `__attrs_post_init__` asserts it for any
  `DataFrameTextFormat` source.

`out_format` is `ParquetOutFormat`. The pipe list may be empty, or carry a column-select pipe
if it is worth dropping `neg_log_10_p_value` / `variant_id` before the index reads it.

### Step 4 — `ConstructCsfVariantIndexTask` + tests

Templated off the aptamer's variants, intersected with the HapMap3 membership list,
allele-aware on the unordered `{EA, NEA}` set; the template's orientation wins.

Output columns: `CHR, POS, EA, NEA, rsID, EAF, is_strand_ambiguous`. **No `POS_HG19`.**
Sorted by the full key `(CHR, POS, EA, NEA)` — row order *is* the alignment slot, so the sort
must be fully deterministic. Sort on all four columns even though `(CHR, POS)` is unique in
today's data; see [File shape](#file-shape).

Expected: **1,025,155 rows**. Assert the count is within a tolerance of that and fail loudly
otherwise; a silent shrink here corrupts every downstream file.

### Step 5 — `BuildSlimCsfAptamerParquetTask` + tests

Simpler than the PPP analogue: no tar, no per-chromosome file loop, no Synapse.

1. Read the index (`CHR, POS, EA, NEA`).
2. `wf.download_from_url` the aptamer's `.tsv.gz` into scratch, md5-verified.
3. Read only `chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, n`.
4. Align to the index, orienting `beta` to the index effect allele; join misses become
   `float('nan')`, **not** polars null (zero-copy numpy conversion downstream).
5. Write `BETA, SE, N` as `float32`, one row group per chromosome, zstd level 15,
   byte-stream-split on all three.

Peak memory is the single-file read, ~0.6–1 GB. That caps useful parallelism; if it proves
too tight, the harmonised `.h.tsv.gz` + `.tbi` variant allows true per-chromosome reads at
the cost of a larger download (291 MB) and a tabix dependency.

Build one aptamer end-to-end and check the output is ~6.8 MB.

### Step 6 — Generator and terminal assets

`generate_csf_slim_aptamer_tasks(index_task, index_name, manifest_path)`, plus the terminal
asset module and the `DEFAULT_RUNNER` entry point. Asset ids follow
`csf_slim_{index_name}_{entrez_gene_symbol}_{seq_id}`; `project` is the gene symbol,
matching PPP's layout.

### Step 7 — Full build

~1.42 TB of transient downloads and ~47.5 GB of output. Downloads dominate: at 50 MB/s
that is ~8 hours of transfer, with ~2 CPU-hours of compression overlapping it.

Route the output subtree to `/mnt/d` via `path_remap` in `default_runner_config.yaml` before
starting — it is exactly the profile that remapping targets (large, many-file, rarely read).

---

## Tests

Following the repo convention of Task-level tests over helper-level ones, and no assertions
on message wording:

| Test | Covers |
|---|---|
| `test_construct_csf_variant_index_task` | Tiny synthetic template + membership; asserts row order, allele orientation adopted from the template, allele-aware matching of a swapped-orientation reference row |
| `test_build_slim_aptamer_parquet_task` | Alignment on a synthetic index: beta sign flip when the effect allele is swapped, NaN for an index variant the aptamer lacks, `N` carried through per variant |
| `test_gwas_catalog_url` | Bucket arithmetic at boundaries (`GCST90421000`, `GCST90421001`, `GCST90422000`) |
| Manifest resolver | The four rules on a small fixture, including one Casein kinase II override |

The `WF` download is injected, per the no-monkeypatching convention — tests supply a fake.

---

## Deferred decisions

**How LDSC consumes `N`.** `constant_sample_size()` in
`mecfs_bio/build_system/task/ppp_ldsc/ppp_protein_heritability_task.py:113` asserts a single
distinct N across context variants. CSF violates that on every aptamer. The options are a
scalar (median) N per aptamer or per-SNP N threaded through `batched_ldsc_h2` /
`batched_ldsc_rg`. Storage no longer constrains this — the data is on disk either way — so it
is a pure correctness question, to be settled when the LDSC tasks are written.

When that happens, the CSF LDSC tasks should be **new task shells** reusing the existing
helper modules (`batched_ldsc_h2`, `batched_ldsc_rg`, `ppp_ldsc_context`, `trait_alignment`),
which import nothing from the PPP task layer and need no changes. Only the two PPP task
shells are concretely typed on `BuildSlimProteinParquetTask`, and they are left alone.

Cis-masking will also need gene coordinates from something other than Sun et al. ST3 — the
Ensembl/MAGMA gene-location assets already in the repo, keyed on UniProt or Entrez symbol.

---

## Risks

| Risk | Mitigation |
|---|---|
| The 47.5 GB estimate assumes all aptamers share one variant set (plausible — one joint PLINK2 run — but unverified). NaN-heavy files would compress *better*, so the estimate is an upper bound on the mean | Check the row counts of the first few built files |
| Catalog re-curation changes a trait string and breaks the resolver | The bijection assertion fails loudly at manifest-regeneration time, not silently at build time |
| 7,008 sequential downloads from EBI FTP; rate limiting or transient 5xx | Use `call_with_retries` (`mecfs_bio/util/retry`), as the OSF task does |
| Aborted mid-build after hours of downloading | The build system is incremental; completed aptamers are not re-fetched |
