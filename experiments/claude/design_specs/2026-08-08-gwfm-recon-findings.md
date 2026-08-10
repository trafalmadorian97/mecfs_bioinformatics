# GWFM Reconnaissance Findings — GCTB / gctbhub

Date: 2026-08-09
Status: Recorded during Task 1 of the SBayesRC GWFM remote-execution plan
(`.superpowers/sdd/2026-08-08-sbayesrc-gwfm-remote-execution/task-1-brief.md`). Feeds
`mecfs_bio/build_system/task/sbayesrc/gctb_gwfm_constants.py`.

All investigation was read-only: `curl -sI` for `Content-Length` (HEAD), directory
listings for filenames, and a single download of the ~2.3 MB binary zip to compute its
sha256. The ~200 GB LD matrix and ~125 GB eigen archive were never downloaded, only
HEAD-probed / listed.

## 1. GWFM reference files + directory URLs

### LD matrix — `.../resources/GWFM/LD/Imputed13M/`

Directory listing:
<https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/>

| File | Size (bytes) | Size | Last-Modified |
|---|---|---|---|
| `ukbEUR_13M_FullLDM.zip` | 206,566,549,726 | 206.57 GB / 192.38 GiB | 2025-08-09 |
| `ref_b37_1588blocks.pos` | 40,058 | 39.1 KiB | 2026-01-13 |
| `eigen/blk1588_eigen.tar.gz` | 125,585,083,484 (from directory listing) | 125.59 GB / 116.96 GiB | 2026-05-15 |

Source: HEAD requests confirmed `ukbEUR_13M_FullLDM.zip` (`Content-Length:
206566549726`) and `ref_b37_1588blocks.pos` (`Content-Length: 40058`) are directly
downloadable, and a `Range: bytes=0-0` request against the zip succeeded
(`Content-Range: bytes 0-0/206566549726`), confirming byte-range support (useful for a
resumable/streaming stager in a later task).

**Anomaly:** `eigen/blk1588_eigen.tar.gz` returns **HTTP 403 Forbidden on both HEAD and
GET** (including `Range` requests), regardless of User-Agent or Referer header — unlike
the sibling `ukbEUR_13M_FullLDM.zip`, which is fully accessible. Its size above is
therefore sourced from the nginx/Apache autoindex directory listing (which the server
itself populates from the file's stat info), not from a successful HEAD. This file could
not be inspected at all. Recommend re-probing before any later task that would download
it.

### Annotation — two candidate directories, only one is correct for GWFM

- `.../resources/v2.0/Annotation/` → `annot_baseline2.2.zip`, 307,548,449 bytes. This is
  the **HapMap3/7M-scale** annotation (same one referenced in the codebase's existing
  `baseline_2_2_annotation` asset per the design doc) — **not** matched to the 13M GWFM
  reference.
  <https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/Annotation/>
- `.../resources/GWFM/Annotation/` → `annot_baseline2.2_13M.zip`, 557,227,563 bytes
  (531.4 MiB), Last-Modified 2026-05-08. This is the **13M-scale** annotation, correctly
  matched to the GWFM LD reference (confirmed by the software page's "13M SNP
  annotations" link pointing at this exact directory).
  <https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/Annotation/>

**Pinned choice: `annot_baseline2.2_13M.zip`.**

### Gene map

`gene_map_hg38_hg19.txt`, 5,155,268 bytes (4.92 MiB), Last-Modified 2026-03-10, linked
directly from the tutorial text ("Gene region annotations for genome build hg19 and hg38
are available here"):
<https://gctbhub.cloud.edu.au/software/gctb/download/gene_map_hg38_hg19.txt>

### Why the precomputed eigen archive is excluded from the pinned bundle

The GWFM tutorial's own worked example (extracted from
<https://gctbhub.cloud.edu.au/software/gctb/> and cross-checked against the `README`
shipped inside `gctb_2.5.5_Linux.zip`) never references `eigen/blk1588_eigen.tar.gz`.
Step 1 always derives a **trait-matched** eigen-decomposition from the **raw** LD matrix:

```
gctb --ldm ukbEUR_13M_FullLDM --gwas-summary test.ma --make-ldm-eigen --thread 32 --out matched_ldm
```

i.e. `--make-ldm-eigen` takes `--ldm` (the raw blockwise full LD matrix folder) plus the
trait's own `--gwas-summary`, and produces `matched_ldm` locally on the compute host —
it is not a downloaded artifact. Step 2 (`--gwfm RC ...`) then consumes that per-trait
`matched_ldm` output via `--ldm-eigen`, not the gctbhub precomputed eigen archive.

The precomputed `eigen/blk1588_eigen.tar.gz` download is instead described on the
software page under a generic, non-GWFM-specific heading: "Eigen-decomposition data of
LD matrices ... for SBayesRC and SBayesR with the low-rank model" — plausibly meant for
the standard (non-genome-wide) SBayesRC/SBayesR workflow, which can skip per-trait
`--make-ldm-eigen` when the trait's SNP set already matches the precomputed one. It sits
in the same `GWFM/LD/Imputed13M/` directory only because it shares the same underlying
13M-SNP reference panel, not because the documented GWFM CLI path uses it.

Combined with the 403 Forbidden on this file (so it can't currently be verified even if
we wanted to use it) and the disk-budget arithmetic in §5 (downloading both the raw
206.57 GB zip *and* the 125.59 GB eigen archive would strain the ~500 GB disk estimate),
**`GWFM_REFERENCE_BUNDLE` deliberately omits `eigen/blk1588_eigen.tar.gz`.** This is a
recorded decision, not a shipped `USES_PRECOMPUTED_EIGEN` boolean — see §3 below and the
ambiguity-resolution note in the task-1 brief.

## 2. Binary

URL: <https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.5.5_Linux.zip>
(`Content-Length: 2294269`, 2.19 MiB, Last-Modified 2026-03-10 — note this Last-Modified
predates the "Last updated: 12 Dec, 2025" banner printed by the binary itself; the
gctbhub file server's Last-Modified header does not reliably track the tool's own
internal release date).

Downloaded once and hashed:

```
sha256sum gctb.zip
ccc2752e1bc0d4a210bd9592d07c44468891677f8b2d5fc0a37aa2119c24c61f  gctb.zip
```

Archive contents: `gctb_2.5.5_Linux/{gctb, LICENSE, README, .DS_Store}` plus
`__MACOSX/` AppleDouble junk. `gctb` is a statically linked 64-bit Linux ELF executable
(`file` confirms "statically linked"), which matters for the minimal `debian-slim`
container image planned in the design doc — no shared-library dependencies to install
beyond what's needed for OpenMP threading (`libgomp1`, used at runtime by `--thread`).

License: `LICENSE` inside the zip is the **MIT License** (Copyright (c) 2017 Jian
Zeng) — confirms the design doc's "already confirmed MIT" note and supports the planned
self-built public image without redistribution concerns.

## 3. Exact GWFM CLI + precomputed-eigen question

The three commands from the starting facts were verified **exactly** three independent
ways: (a) the software page's rendered tutorial text, (b) the plain-text `README` shipped
inside the binary zip, and (c) `strings` on the compiled binary itself, which contains
every flag name used below (`--ldm`, `--gwas-summary`, `--make-ldm-eigen`, `--thread`,
`--out`, `--gwfm`, `--ldm-eigen`, `--annot`, `--gene-map`, `--cs`, `--pwld-file`, `--pip`,
`--pep`, `--mcmc-samples`, `--flank`, `--genome-build`, `--block`, `--get-pwld`).

```
# Step 1 — match + eigen-decompose the LD reference to the trait's SNPs
gctb --ldm ukbEUR_13M_FullLDM --gwas-summary test.ma --make-ldm-eigen --thread 32 --out matched_ldm

# Step 2 — genome-wide fine-mapping (the ~13 h / ~150 GB heavy step)
gctb --gwfm RC --ldm-eigen matched_ldm --gwas-summary test.ma --annot annot.txt --gene-map gene_map.txt --thread 32 --out test

# Step 3 — (re)calculate credible sets (cheap)
gctb --cs --pwld-file ldm/rsq0.5.pwld --pip 0.9 --pep 0.7 --gene-map gene_map.txt --flank 5000 --genome-build hg19 --mcmc-samples test --out test
```

Source: <https://gctbhub.cloud.edu.au/software/gctb/> (rendered "Genome-wide
fine-mapping (GWFM)" tutorial section) and the `README` inside
`gctb_2.5.5_Linux.zip`.

Note: `GCTB_CS_TEMPLATE` in the constants module omits `--flank 5000 --genome-build
hg19` because both match the tool's documented defaults (flank default 5000, genome
build default hg19); they're not required flags for the pinned use case (hg19/GRCh37
reference throughout this plan). If a later task needs hg38, add `--genome-build` to the
template's substitutable fields rather than hardcoding it.

**Does `--gwfm RC` consume the precomputed `eigen/` directory directly, skipping Step
1?** No evidence found that it does — see §1's "Why the precomputed eigen archive is
excluded" analysis above for the full reasoning. This is recorded as a finding only, per
the task-1 brief's ambiguity resolution #2: **Task 10 (the GWFM task) should always emit
Step 1 (`--make-ldm-eigen`)**; it should not attempt to shortcut it using a shared
precomputed eigen bundle. No `USES_PRECOMPUTED_EIGEN` constant was added.

## 4. Checkpoint / resume support

**No evidence of MCMC checkpoint/resume support in the `gctb` binary itself.**

- The GWFM tutorial text and binary `README` contain no mention of resume, checkpoint,
  or `.rds` files.
- `strings` on the `gctb` executable turned up no resume/checkpoint-related flags or
  messages — only an unrelated C++ runtime symbol (`_Unwind_Resume`) and one MCMC
  diagnostic string, "Restarting MCMC with a more robust parameterisation for SBayes"
  (a convergence-recovery restart *within* a single run, not a resume-after-interruption
  feature).
- `gctb --help` is not a recognized flag ("Error: invalid option \"--help\""); running
  with no arguments just prints the banner and "Did you forget to give the input
  parameters?" — there is no built-in flag listing to double-check against, so absence
  of resume flags in `strings` is the best available evidence, not a guaranteed
  negative.

This is a **different tool** from the "SBayesRC silently resumes from
`<outPrefix>.rds`" behavior noted elsewhere in project memory
(`project_ppp_sbayesrc_feasibility` / polypwas work) — that note is about the
**SBayesRC R package** used for polygenic-prediction (`zhiliz/sbayesrc` Docker image,
`Rscript` entry point), not the pure C++ `gctb` binary used here for GWFM. No
`GCTB_SUPPORTS_RESUME` constant was added; per the brief, this stays out of scope unless
revisited later (it currently gates the design doc's deferred spot-instance decision
toward "do not use spot" since there's no confirmed resumption path for `gctb --gwfm`).

## 5. VM disk sizing → `DEFAULT_DISK_GB`

Known artifact sizes relevant to the disk budget:

| Item | Size |
|---|---|
| `ukbEUR_13M_FullLDM.zip` (download) | 206.57 GB |
| `ukbEUR_13M_FullLDM` unzipped | unknown; not measured (never downloaded) |
| `annot_baseline2.2_13M.zip` (+ unzipped) | 0.56 GB + similar unzipped |
| `gene_map_hg38_hg19.txt` | 5.2 MB |
| Step-1 `matched_ldm` eigen output | unknown; expected well below the 206.57 GB raw matrix since it's a low-rank, block-restricted decomposition, but not measured |
| Step-2/3 result files (`.snpRes`, `.mcmcsamples`, `.parRes`, etc., ~13M SNPs) | unknown; not measured |

Because the raw LD matrix is stored as dense per-block correlation data, zip compression
on it is expected to be modest (floats don't compress well) — plausibly 1.0–1.3×
unzipped-to-zipped, i.e. **~210–270 GB unzipped**. Simultaneously holding the 206.57 GB
zip *and* its unzipped form would already consume ~420–480 GB; the design assumes the
staging/setup step deletes (or streams-and-discards) the zip after extraction rather than
keeping both. Add the annotation (~1 GB unzipped), gene map (negligible), the Step-1
eigen output and Step-2/3 result files (unmeasured, but expected in the tens-of-GB range
based on 13M-SNP-scale per-SNP result tables), plus OS/Docker image overhead (~10–20
GB).

Under the "delete zip after unzip" assumption, ~500 GB has comfortable headroom (~230 GB
peak for the raw matrix + double-digit GB for everything else); under a
"keep zip and unzipped copy simultaneously" assumption it would not. **This is an
estimate, not a measurement** — the brief's ambiguity resolution already fixes
`DEFAULT_DISK_GB = 500`, which this reconnaissance supports but a later staging/setup
task should verify empirically (actual unzip ratio, actual eigen/result sizes) and bump
if needed.

## 6. Zip inner structure (probed 2026-08-09 for Task 10, via HTTP Range on each zip's central directory — no full download)

Both zips are single-root archives. Read by range-fetching the tail (End-Of-Central-Directory + central directory) of each file; see experiments/claude/logs/gwfm_zip_central_dir_probe.log.

- `annot_baseline2.2_13M.zip` (1 entry) unzips to a single file: **`annot_baseline2.2_13M.txt`**. So GWFM step 2's `--annot` path is `annot_baseline2.2_13M.txt` (NOT the tutorial's generic "annot.txt").
- `ukbEUR_13M_FullLDM.zip` (1592 entries, zip64) unzips to a single folder **`ldm13M/`** (NOT "ukbEUR_13M_FullLDM" as the tutorial example's `--ldm` arg loosely suggested). Its non-block members are `ldm13M/snp.info`, `ldm13M/ldm.info`, `ldm13M/rsq0.5.pwld`, plus 1588 `ldm13M/block<N>.ldm.bin`. So:
  - GWFM step 1 `--ldm` (raw blockwise LD folder) = `ldm13M`
  - GWFM step 3 `--pwld-file` = `ldm13M/rsq0.5.pwld` (it ships inside the LD folder; no separate `--get-pwld` step needed)
- `ref_b37_1588blocks.pos` and `gene_map_hg38_hg19.txt` are staged as-is (no unzip); `gene_map_hg38_hg19.txt` is the `--gene-map` arg. `ref_b37_1588blocks.pos` is NOT referenced by any of the three CLI templates — staged for completeness / potential future `--block` use, flagged.

These inner names are now pinned as constants consumed by Task 10's command construction; a toy Task-12 reference must reproduce the same names (a `ukbEUR_13M_FullLDM.zip` unzipping to `ldm13M/` with `snp.info`/`ldm.info`/`rsq0.5.pwld` + toy blocks, and an `annot_baseline2.2_13M.txt`).

## Open items / concerns for later tasks

1. `eigen/blk1588_eigen.tar.gz` is 403-Forbidden on every request tried (HEAD, GET,
   Range, with/without UA and Referer headers). Not blocking for this task (it's
   excluded from the bundle — see §1), but worth a sanity re-check before any future
   task considers using it.
2. `DEFAULT_DISK_GB = 500` is derived from an estimate (raw zip + unzip ratio guess), not
   a measured unzip size, matched-eigen size, or result-file size. Flag for verification
   during the first real staging/GWFM run.
3. `GWFM_REFERENCE_VERSION = "Imputed13M/v1"` is the value fixed by the task-1 brief; it
   is not a version string gctbhub itself publishes (the site has no explicit
   "Imputed13M v1/v2" versioning visible in directory listings) — treat it as *our*
   pinned label for this snapshot of the reference files, to be bumped if/when gctbhub
   republishes any of the pinned files.
