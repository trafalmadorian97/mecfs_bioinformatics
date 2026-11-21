# deCODE SomaScan proteomics database and coloc window reader — design

**Date:** 2026-11-21
**Status:** draft, awaiting user review
**Branch:** `decode-somascan-database`

## Intent

The repo has proteomics "databases" for UKB-PPP (`mecfs_bio/assets/gwas/ukbb_ppp/ppp_database`)
and the Western et al. CSF pQTL study (`mecfs_bio/assets/gwas/csf_pqtl/csf_database`). Each is a
shared HapMap3 variant index plus one slim beta/se(/N) parquet per protein, aligned row for row to
the index. They feed genetic correlation (PPP LDSC pipeline) and the hit table used to screen
proteins for colocalization. Coloc itself (`protein_coloc/trait_protein_coloc_task.py`) reads
dense summary statistics for a few locus windows per protein through a source-specific
`ProteinWindowReader` (`PppWindowReader`, `CsfWindowReader`).

We want the same for the **Icelandic SomaScan plasma proteomics GWAS from deCODE** (Eldjarn et
al. 2023, "Large-scale plasma proteomics comparisons through genetics and disease
associations"; the SomaScan arm is Iceland, n ≈ 30–36k per variant). A large Icelandic plasma
SomaScan study complements UKB-PPP (Olink, UK) and Western (CSF, SomaScan 7k).

Scope is the **SomaScan data only** (user decision). The bucket also holds deCODE's UK Biobank
Olink reanalyses (British/Irish, African, South Asian); those are out of scope.

## What the survey established

Survey scripts and logs: `experiments/claude/decode_proteomics_survey/` (logs not committed).

**Access.** An S3-compatible bucket, not a temporary URL. deCODE emails a bucket name
(`largescaleplasma-2023`), endpoint (`s3-ext.decode.is:443`, HTTPS) and an access key pair that
expires (current one: 18 Dec 2026). Expiry is not a deadline: the user can request new keys
whenever needed, so the design just needs credentials that are easy to rotate. deCODE's
instructions (`aws configure --profile decode-download`, region `us-east-1`) recommend a named
AWS profile. The server is NetApp StorageGRID 12.1.

**Layout.**

| Prefix | Files | Size | Meaning |
|---|---|---|---|
| `final_somascan_smp/` | 5,284 | 5.13 TB | Iceland SomaScan, SMP-normalized (`Proteomics_SMP_PC0_<seq>_<seq2>_<GENE>_<name>_10032022.txt.gz`) |
| `final_somascan_raw/` | 5,284 | 5.13 TB | Iceland SomaScan, not SMP-normalized (`Proteomics_PC0_..._07082019.txt.gz`) |
| `variant_info_iceland.txt.gz` | 1 | 442 MB | `Chrom Pos A1 ImpFreqA1 ImpInfo` |
| `proteomics2023README_20240517.txt` | 1 | 1 KB | column definitions |

Both groups list the same 5,284 SomaScan seq ids, one file per seq id.

**File format.** Each file (~970 MB) is one plain gzip stream, not BGZF and with no tabix index,
holding tab-separated text (3.32 GB uncompressed, 33,449,999 variants, GRCh38):
`Chrom Pos Name rsids effectAllele otherAllele Beta Pval minus_log10_pval SE N ImpMAF`.
- `Chrom` is `chr1`…`chrX`. Rows are sorted by position within a chromosome, but **chromosomes
  are in text order** (chr1, chr10, …, chr19, chr2, chr20, …, chr9, chrX).
- `N` varies per variant (~30.4k–35.7k).
- `ImpMAF` is the minor allele frequency; there is no effect-allele frequency column.
- `rsids` uses `.` for missing.
- A multi-allelic test is written with `!` in `otherAllele` (meaning "not effectAllele"). None
  appeared in the profiled file, but the README documents them.
- Each object carries deCODE's md5 in its metadata (`x-amz-meta-md5chksum`, base64).
- HTTP Range requests work (206).

**Shared row order.** Every SomaScan file appears to list the same variants in the same order:
- AP4M1 and CRYBB2 (both SMP) have identical full-length row counts.
- All 1,197,450 HapMap3 rows picked by row number from AP4M1 match CRYBB2's variant names.
- The first ~864k variant names are identical in 7 files across both groups.

This is not guaranteed by deCODE, so the per-protein task checks it (see below).

**Throughput.** deCODE caps each HTTP connection at ~0.6 MB/s (≈5 Mbit/s, constant regardless of
conditions), but the cap is per request, so throughput grows with connections up to an overall
ceiling of ~25–30 MB/s:

| Connections | Total |
|---|---|
| 1 | 0.7 MB/s |
| 8 | 4.9 MB/s |
| 32 | 15.5 MB/s |
| 64 | 21.8 MB/s |
| 128 | 26.8 MB/s |

The ceiling is not this machine's link (≥ 46 MB/s to another host). One file downloaded with 128
parallel range requests took 38–43 s (22–27 MB/s), and its md5 matched.

**Mockup** (`mockup_single_protein.py`, template + protein modes): the full per-protein build
takes **53.8 s** (43 s download, 2 s md5, 8 s processing), with peak memory 4.2 GB. Remote window
reads take 1.4–1.7 s and fetch ~6 MB each, and both matched a brute-force scan of the file.

## Goals

1. **Credentials** in a user-local, uncommitted place that is easy to rotate. Tasks never contain
   keys or URLs that expire.
2. A **deCODE SomaScan database per group** (SMP and non-normalized), each built the same way:
   - a shared HapMap3 variant index;
   - one slim per-protein parquet (BETA, SE, N; float32; index row order), usable by the same
     consumers as the PPP/CSF slim files: hit table and, later, LDSC;
   - a per-protein **gzip seek index** and **position map**, so coloc can read dense windows from
     deCODE's remote files without downloading them.
3. A **hit table** asset per group (p < 5e-8, as for PPP/CSF).
4. A **`DecodeWindowReader`** implementing `ProteinWindowReader`, plus the coloc wiring (asset
   generator group, gene-coordinate and UniProt-accession readers) so `TraitProteinColocTask`
   can run against deCODE like PPP and CSF.
5. The serial build system as it is, with no concurrent-task support, and **resumable** builds:
   one task per protein, so an interruption costs at most one protein.

## Non-goals

- deCODE's UK Biobank Olink files.
- Running coloc or LDSC for specific traits. That is follow-up work once the database exists;
  this design only makes it possible. *(Proposal; see open decisions.)*
- Keeping any full dense summary-statistics file on disk.
- Concurrent task execution in the build system.
- Any text under `docs/`.

## Anti-goals (a failure even if it runs)

- Committing credentials or any expiring URL, or putting either in a Task's fields.
- A design that only works while one set of keys is valid. Once built, the database (index,
  slim, seek index, position map) is local and usable without credentials. Only coloc window
  reads need valid keys.
- Silently assuming the shared row order. A protein whose rows differ must fail its task loudly.
- A long-running batch task whose interruption restarts the whole download.

## Constraints

- Repo conventions (CLAUDE.md and the user's memory notes), especially:
  - attrs frozen classes;
  - Literal/NewType ids;
  - external services via `WF`, injected and faked in tests (no monkeypatching);
  - `call_with_retries` as the one backoff loop;
  - polars;
  - specifics (bucket, prefixes, profile name) under `assets/`, classes in `build_system/`;
  - new dependencies in the pixi analysis feature;
  - Task-level unit tests;
  - imports at top of file.
- Build machine: 16 cores, 15 GB RAM (WSL). The 16 GB CI runner is an intentional limit.
- Network: ~25–30 MB/s overall, at most ~0.6 MB/s per connection.

## Approaches considered

**Getting around the slow download.**

| Option | Verdict |
|---|---|
| aria2 (repo's `AriaDownloader`) | Capped at 16 connections per server, so ~8 MB/s per file. It needs ~8 files downloading at once, which the serial build can't do. |
| One batch task with its own concurrency | Fast, but an interruption days in restarts it. Rejected (anti-goal). |
| One task per chunk of proteins | Workable middle ground, but unnecessary given the next row. |
| **One task per protein, with a downloader that splits one file across ~128 range requests** | **Chosen.** s5cmd (19.5 MB/s), rclone (21.4) and a plain Python thread pool over `urllib` (24.2) all exceed aria2's cap. Python needs no new tool and is easy to fake in tests. |
| Cloud VM near Iceland | Not needed. It only helps if the ceiling is per client IP. |

**Coloc window access.**

| Option | Verdict |
|---|---|
| Tabix | Impossible: tabix needs BGZF, and deCODE's files are plain gzip, which we can't change. |
| Local dense mirror (re-compressed BGZF) | About 5 TB per group. Rejected. |
| Stream from the start of the file to the window | About 500 MB per protein per window (~100 s with 8 connections). Rejected. |
| **zran-style gzip seek index built during the per-protein pass, plus a position map** | **Chosen.** About 1 MB per protein; reads fetch ~6 MB per window in ~1.5 s. |

**Building the seek index.**

| Option | Verdict |
|---|---|
| `indexed_gzip` sequential read | 27 s per file. |
| **`rapidgzip`** (8 threads, `chunk_size=16 MiB`), exporting `IndexFormat.INDEXED_GZIP` | **Chosen.** 5 s per file, a smaller index (0.92 MB zstd, 232 points, sparse windows), and `indexed_gzip` loads it for reads. Verified by random seek/read against plain zlib. |

**Picking the HapMap3 rows.**

| Option | Verdict |
|---|---|
| Full parse plus a join on chromosome, position and alleles | 9 s per file (join 5–7 s); integer-key and position-first variants were no faster. |
| **Parse only the shared index's rows, located by row number with a numpy newline scan** | **Chosen.** 3 s per file, and the result was identical. It relies on the shared row order, which each task checks. |

## Design

### 1. Credentials and the S3 capability (WF)

- **Credentials** live in a named AWS profile, `decode-download`, in `~/.aws/credentials` and
  `~/.aws/config`, exactly as deCODE's instructions set it up. The config gets
  `region = us-east-1`, and `endpoint_url = https://s3-ext.decode.is` is optional; we pass the
  endpoint explicitly. Rotating keys means rerunning `aws configure --profile decode-download`.
  The profile name, endpoint and bucket are constants under `assets/` (they are deCODE
  specifics).
- **New WF capability** `s3_access: S3Access` (Protocol in `build_system/wf/`), defaulted in
  `make_wf`. It names a profile, endpoint and bucket and provides:
  - `head(key) -> S3ObjectInfo` (size, deCODE md5 from `md5chksum`);
  - `presigned_url(key, expires_s) -> str`.

  The production implementation uses `boto3` with SigV4 and path-style addressing. Tasks get
  presigned URLs at execution time only.
- **Failure messages.** A missing profile, or a 403 from expired keys, fails fast with a message
  naming the profile and the `aws configure` command.
- **New dependencies:** `boto3`, `rapidgzip` (0.16.2), `indexed_gzip` (1.10.3) and `zstandard`
  (0.25.0), all in the pixi analysis feature. All four install from conda-forge for Python 3.13
  (checked). The system `aws` CLI is not used by production code. Python 3.13 has no stdlib
  zstd, hence `zstandard`.

### 2. Parallel range downloader

`download_with_parallel_ranges(url, size, out_path, n_connections=128, chunk_bytes=4 MiB, ...)`
in `util/download/`:
- splits the object into chunks and fetches them with a thread pool of HTTP Range GETs, each
  wrapped in `call_with_retries` (retrying connection errors, timeouts and 5xx);
- `pwrite`s each chunk at its offset into a preallocated file;
- asserts 206 responses and exact lengths.

The md5 is checked after the download against `S3ObjectInfo`. It is exposed on `WF` (e.g.
`wf.download_url_parallel(...)`) so tests inject a fake that copies a local file. The defaults
come from the survey and are tunable. Smaller chunks (2–4 MiB) are expected to beat 8 MiB
slightly, so the builder may measure and pick.

### 3. Manifests (committed)

A one-off regeneration script (like `regenerate_ppp_manifest.py`) lists both prefixes, sends a
HEAD to each object (for size and md5), parses file names, and writes one committed CSV per
group under `mecfs_bio/assets/reference_data/decode_somascan_sumstats/`. Columns:
- `seq_id`, in dash form (`10000-28`, from `Proteomics_*_PC0_10000_28_…`);
- `gene_symbol` from the file name;
- `target_name`, the rest of the file name;
- `key`, `size`, `md5`;
- `uniprot` and `entrez_gene_id`, from the annotation below.

**Protein annotation** (UniProt and Entrez per seq id): the CSF manifest
(`csf_aptamer_manifest.csv`, SomaScan 7k) already covers **4,813 of the 5,284** deCODE seq ids.
The rest need deCODE's own aptamer table (Ferkingstad et al. 2021 / Eldjarn et al. 2023
supplementary tables). The builder reads those papers' supplements (per repo practice: read
the PDF, don't work from memory), uses deCODE's table as the primary source, and records
disagreements with the CSF manifest. Seq ids left without UniProt get no accession rows and an
"unknown" cis status in coloc, as CSF already handles.

### 4. Variant info asset (effect-allele frequency)

The sumstats only carry `ImpMAF`. `variant_info_iceland.txt.gz` (442 MB) has `ImpFreqA1`.
- **Task:** download it once (parallel range downloader + md5) and convert it to parquet
  (`Chrom` code, `Pos`, `A1`, `ImpFreqA1`, `ImpInfo`). This is cheap and needs credentials only
  once.
- **Index construction:** join it in on chromosome, position and effect allele. Assert that
  `A1 == effectAllele` for every index row; only one row was spot-checked so far. If that holds,
  the index has a true `EAF`. If it doesn't, the builder stops and reports: we would fall back to
  MAF and name the column honestly.

### 5. Shared variant index (one per group)

`ConstructDecodeVariantIndexTask`, modelled on `ConstructCsfVariantIndexTask`.
- **Template:** one template protein file (the first manifest row of the group, so the choice
  is deterministic), downloaded with the parallel range downloader, decompressed with
  `rapidgzip`, fully parsed, and joined to `HAPMAP_3_MEMBERSHIP_LIST` on chromosome, position
  and the unordered allele pair. The deCODE orientation is adopted. The template protein's own
  per-protein task still runs normally; it is checked against the index like every other
  protein.
- **Columns:**
  - the usual index columns: CHR, POS (hg38), EA, NEA, rsID, EAF from §4;
  - `SOURCE_ROW` (UInt64), the variant's 0-based data row in every deCODE file of the group;
  - `NAME`, deCODE's `chr:pos:ref:alt`, used for checking.
- **Parquet metadata:** the group's `total_rows` (33,449,999 for SMP).
- **Rows excluded:** `!` (multi-allelic) rows.
- **Output:** 1,197,450 rows, 11.3 MB in the mockup; about 60 s to build.

### 6. Per-protein task: `BuildDecodeSomascanProteinTask`

Inputs: the group's index task, the manifest row, the S3 key and md5. In `execute`:
1. HEAD the object, check its size and md5 against the manifest, get a presigned URL, and
   download with parallel ranges into `scratch_dir`. Check the md5.
2. **One `rapidgzip` pass** (`parallelization=8`, `chunk_size=16 MiB`, reading 16 MiB blocks):
   - **Line boundaries:** a numpy newline scan per block gives row numbers and offsets.
   - **Position map:** a row for the first line starting at or after each 1 MiB of
     uncompressed text: `chrom` Int8 (X=23), `pos` Int64, `line_offset` UInt64, `row_number`
     UInt64. That's about 3,300 rows.
   - **HapMap3 rows:** slice out only the `SOURCE_ROW` lines and parse them with polars.
   - **Checks:**
     - the header matches the expected columns;
     - the total row count equals the index's `total_rows`;
     - every gathered `Name` equals the index's `NAME`;
     - no gathered `otherAllele` is `!`.
   - After the pass, export the seek index with `export_index(..., IndexFormat.INDEXED_GZIP)`.
3. **Output: a `DirectoryAsset` with three files:**
   - `slim.parquet`: BETA, SE, N as float32 in index row order (zstd plus byte-stream-split,
     per the docs-parquet convention; 7.9 MB in the mockup, with N per variant as for CSF);
   - `seek_index.gzidx.zst`: the GZIDX v1 export, zstd-compressed (0.92 MB);
   - `posmap.parquet`, whose parquet metadata holds the source key, compressed size, md5, the
     header's end offset, and the format versions.
4. The raw download is deleted along with the scratch directory.

Budget per protein: about 54 s wall time, about 8.8 MB on disk. Peak memory was 4.2 GB in the
mockup. The builder should bring it down, for example with fewer `rapidgzip` workers, smaller
blocks and no 16 MiB buffer concatenation, and record the trade-off. A target of 2 GB or less
is my proposal.

**Consumers of the slim file.** `scan_dataframe_asset` requires a `FileAsset`, but
`ProteinHitTableTask` (and the PPP LDSC tasks) read `task` outputs through it. The builder adds
the smallest generic way to read a dataframe stored as a named file inside a `DirectoryAsset`,
for example a read-spec or meta field naming the member file, which `scan_dataframe_asset`
honours. The hit table then works unchanged with deCODE tasks. Copying the slim file into its
own task is rejected because it duplicates ~42 GB.

### 7. Asset generators and asset layout

- **Generator:** `generate_decode_somascan_protein_tasks(group, index_task, index_name)` returns
  a collection of per-protein tasks, one per manifest row. `group` is a
  `Literal["smp", "raw"]`.
- **Asset ids and paths** are namespaced by group and index, e.g.
  `decode_somascan_smp_hapmap_3_<GENE>_<seq_id_token>`.
- **Modules:** `mecfs_bio/assets/gwas/decode_somascan/` gets `raw/` (manifest-driven
  specifics), `decode_somascan_database/hapmap3/` (index, protein files, hit table), and a coloc
  asset-generator group mirroring CSF's.
- Both groups are defined. The user picks which to build first (open decision).

### 8. Hit table

`ProteinHitTableTask.create(...)` over each group's per-protein tasks, with
`protein_id = ProteinId(seq_id)` and `gene = gene_symbol`, at p < 5e-8, like the CSF hit table
(`hapmap3_csf_hit_table.py`).

### 9. `DecodeWindowReader` (coloc)

Implements `ProteinWindowReader`. It is built with:
- the group's per-protein tasks (to find each protein's seek index and position map in the local
  asset store, via `fetch`);
- the S3 key per protein;
- the variant info asset (for EAF).

`fetch_protein(protein_id, windows, scratch_dir, wf)`:
1. **Load** the protein's position map and seek index. Decompress the index into memory;
   `indexed_gzip.import_index` accepts a file object.
2. **Validate:** HEAD the object and check its size and md5 against the position map's metadata.
   The index is only valid for those exact bytes, so a mismatch is an error, not a retry.
3. **For each window**, in file order:
   - find the uncompressed span `[U0, U1)` in the position map. Searches go by file order (text
     order of chromosomes), never by chromosome number. The begin is the last row on the
     chromosome at or before `start` (or the row before that chromosome's first row). The end
     is the first later row past the chromosome's section or with `pos > end`, else EOF.
   - find the compressed span from the seek points: from the point at or before `U0` minus one
     byte, to the first point after `U1` plus a 1 MiB read buffer.
   - **prefetch** the compressed span with ~16 parallel range requests of ~512 KB, since each
     connection gets only ~0.6 MB/s.
   - **decompress:** run `IndexedGzipFile(fileobj=<remote source>, index_file=…,
     auto_build=False)`, then `seek(U0)` and `read(U1 - U0)`. The remote source is a seekable
     file-like object serving the prefetched bytes and range-reading anything else.
   - **parse** and keep only rows inside the window.
4. **Convert to `WINDOW_SCHEMA`:**
   - `chrN` becomes an integer, with X = 23;
   - effectAllele becomes EA and otherAllele becomes NEA;
   - Beta, SE and N map across;
   - EAF is joined from the variant info asset on chromosome, position and EA;
   - `!` rows are dropped.

   Orientation to the reference genome stays downstream, as for every source.
5. **Concurrency:** wrap it in `ParallelWindowFetcher` with `DECODE_MAX_PARALLEL_FETCHES`. The
   builder measures this, given the ~0.6 MB/s per connection and the ~25–30 MB/s ceiling;
   for example 4–8 proteins times 16 connections.

### 10. Coloc wiring

- **A gene-coordinate reader and a UniProt-accession reader for deCODE.** The CSF manifest keys
  gene coordinates by Entrez id against the GRCh38 MAGMA table, with the same conflict rules as
  `csf_gene_coords_reader.py`. The deCODE readers should reuse that logic, generalised if
  needed.
- **An asset-generator group** (like CSF's) producing the three coloc outputs for a trait, plus
  a driver module ready for `run_*` functions per trait. Running traits is follow-up work.

## Performance and storage budget (per SomaScan group)

| | Per protein | 5,284 proteins |
|---|---|---|
| Build time (download-bound) | ~54 s | ~3.3 days |
| Slim parquet | 7.9 MB | ~42 GB |
| Seek index (zstd) | 0.92 MB | ~4.9 GB |
| Position map | 0.03 MB | ~0.2 GB |
| Coloc window read | ~1.5 s, ~6 MB | — |

The build will hit transient errors and expiring presigned URLs. Presign a URL per task (6 h
expiry is ample for one protein); each task is independent and resumable.

## Testing

All unit tests run offline:
- `S3Access` and the parallel downloader are injected fakes that serve a synthetic local gzip
  file.
- The synthetic data is a small file in deCODE's format: several chromosomes in text order, a
  `!` row, multiple blocks per seek point (small `chunk_size` and spacing).
- **Index task:** HapMap3 rows, `SOURCE_ROW`, `total_rows`, EAF join and the A1 assertion.
- **Per-protein task:** output files and values against a brute-force parse. It fails loudly on
  a shuffled row, a wrong total, or an md5 mismatch.
- **`DecodeWindowReader`:** windows at the start and end of the file, spanning a seek point, on
  a short chromosome and on chrX, compared with a brute-force filter. A fake range source counts
  bytes, so the "prefetch covers the reads" property is checked.
- **Directory-asset dataframe read** (§6): a small build-system test.
- **Index compatibility round trip:** a `rapidgzip` export, imported by `indexed_gzip`, then
  seeks and reads at random offsets, compared with plain zlib. This guards the two libraries'
  format compatibility, which the design relies on, against version bumps.

No network in unit tests. The real-bucket check is a documented manual smoke run (one protein
end to end, mirroring the mockup).

## Risks

- **Row order changes in some file.** Caught per task. Recovery would be a per-protein fallback
  to a full parse and join, which isn't built unless it happens.
- **deCODE changes or re-uploads a file.** The md5 and size checks in both the build and the
  window reader catch it. Rebuild that protein.
- **Throttling changes.** The downloader's connection count is a parameter.
- **`rapidgzip` and `indexed_gzip` maturity and format compatibility.** Pin versions. The
  round-trip property (rapidgzip export, indexed_gzip import, seek and read) gets a unit test.
- **Credential expiry mid-build.** Tasks fail fast with a renewal message, and completed
  proteins are kept.
- **Annotation gaps** for the ~471 seq ids missing from the CSF manifest (§3).

## Left to the builder

- Exact defaults: download chunk size and connection count, `rapidgzip` workers and chunk size,
  `DECODE_MAX_PARALLEL_FETCHES`.
- Reducing memory, against the target above.
- The exact form of the directory-member dataframe read (§6), kept minimal and generic.
- Module and class names beyond those given here.

## Open decisions for the user

1. **Which SomaScan group to build first** (SMP or not normalized). Both are defined; building
   one takes about 3.3 days.
2. **Scope:** the proposal stops at a working `DecodeWindowReader` plus coloc wiring, with no
   trait runs and no LDSC. Should a first trait run (e.g. DecodeME) or LDSC h²/rg over the deCODE
   database be included?
