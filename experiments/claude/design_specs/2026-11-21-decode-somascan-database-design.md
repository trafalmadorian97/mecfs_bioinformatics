# deCODE SomaScan proteomics database and coloc window reader — design

**Date:** 2026-11-21
**Status:** draft, under user review
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

Both groups list the same 5,284 SomaScan seq ids, one file per seq id. Eldjarn et al. analysed
only **4,907** of them. Their supplementary table ST29 lists exactly those 4,907 seq ids, and
all 4,907 are in the bucket. The other 377 are presumably the non-human and defective aptamers
the paper excludes. The database covers the paper's 4,907 assays only (§3).

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

## Choice of SomaScan normalization: SMP first

Sources: Eldjarn et al. 2023 main text, Supplementary Note 4, and supplementary tables ST19,
ST20, ST29 and ST36. They were read on 2026-11-21; see the project memory
`project_decode_somascan_smp_vs_pc0`.

**What SMP normalization is.** SomaLogic's manufacturer-recommended last step. Each sample is
rescaled to a reference, using one scaling factor per dilution group. Ferkingstad et al. 2021
skipped it ("non-normalized", files `Proteomics_PC0_*`). Eldjarn et al. ran every analysis both
ways and report SMP in the main text.

| | Non-normalized | SMP |
|---|---|---|
| Individuals | 36,136 | 35,892 |
| Assays with a cis pQTL | 1,889 | **2,120** |
| Median top cis −log10 p, assays with a cis pQTL in both (1,880) | 83 | **131** |
| Median correlation between assays (within platform) | 0.42 | −0.08 |
| Median correlation with matching Olink assays | 0.39 | 0.33 |

The cis rows were computed from ST19/ST20. Of the assays with a cis pQTL, 1,880 have one in both
versions, 240 only under SMP and 9 only non-normalized. SMP gives the stronger signal for 1,685
of the 1,880.

**Decision (user, 2026-11-21): build the SMP group first.** Coloc is driven mostly by cis pQTLs,
and SMP has clearly more cis power. The authors say SMP "tends to increase" power for
associations that are not pleiotropic.

**What SMP costs, carried into the risks below:**
- It removes a shared global factor of protein levels. The authors suggest part of that is
  biological.
- It reshapes trans signals at pleiotropic variants. The ASGR1 del12 variant raises 556 proteins
  and lowers none in non-normalized data, but raises 288 and lowers 331 under SMP.
- 18% of protein–phenotype associations significant in both versions flip direction, mostly for
  BMI, blood counts, inflammation, and kidney and liver traits.

The **non-normalized group** is defined by the same code and built later, as a sensitivity
analysis. It is the natural check for genetic correlation (both versions are affected by the
global factor, in opposite ways) and for any important trans coloc.

## Goals

1. **Credentials** in a user-local, uncommitted place that is easy to rotate. Tasks never contain
   keys or URLs that expire.
2. A **deCODE SomaScan database per group** (SMP and non-normalized), each built the same way.
   The SMP group is built first; the non-normalized one is defined but built later as a
   sensitivity analysis. Each group has:
   - a shared HapMap3 variant index;
   - one slim per-protein parquet (BETA, SE, N; float32; index row order), usable by the same
     consumers as the PPP/CSF slim files: hit table and LDSC heritability;
   - a per-protein **gzip seek index** and **position map**, so coloc can read dense windows from
     deCODE's remote files without downloading them.
3. A **hit table** asset per group (p < 5e-8, as for PPP/CSF).
4. A **`DecodeWindowReader`** implementing `ProteinWindowReader`, plus the coloc wiring (asset
   generator group, gene-coordinate and UniProt-accession readers) so `TraitProteinColocTask`
   can run against deCODE like PPP and CSF.
5. The serial build system as it is, with no concurrent-task support, and **resumable** builds:
   one task per protein, so an interruption costs at most one protein.
6. **LDSC heritability** (all-variants and cis-excluded) for every protein of the SMP group, as
   the PPP pipeline produces for UKB-PPP (§11).
7. An **exploratory DecodeME × deCODE coloc run** at trait lead p < 1e-5, with selection
   budgets set so that **no trait-locus × protein candidate pair is dropped for budget** (§12).

## Non-goals

- deCODE's UK Biobank Olink files.
- Genetic correlation (rg) between traits and deCODE proteins, LCV, and coloc runs for traits
  other than the exploratory DecodeME run. These are follow-ups.
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
  - specifics (endpoint, bucket, prefixes, profile name) under `assets/`, classes in
    `build_system/`;
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

### 1. Credentials and S3 access

**Split (user decision 2026-11-21, "option A").** Tasks name *what* they fetch: endpoint,
bucket, key and the name of the credentials profile. `WF` supplies *how* to reach any S3
endpoint. This mirrors Synapse, where tasks hold a `synid` and `WF`'s
`SharedClientSynapseDownloader` logs in once for all of them. `WF` holds no bucket or endpoint,
so one run can read any number of buckets.

- **Credentials** live in a named AWS profile, `decode-download`, in `~/.aws/credentials` and
  `~/.aws/config`, exactly as deCODE's instructions set it up (region `us-east-1`). Rotating keys
  means rerunning `aws configure --profile decode-download`. The profile name is a public label,
  not a secret, so it can sit in task fields.
- **Value types** (`build_system/`, attrs frozen):

  ```python
  @frozen
  class S3Endpoint:
      url: str                  # "https://s3-ext.decode.is"
      region: str               # "us-east-1"
      credentials_profile: str  # "decode-download": AWS profile name, not a secret

  @frozen
  class S3ObjectRef:
      endpoint: S3Endpoint
      bucket: str
      key: str
  ```

  The deCODE values (`DECODE_S3_ENDPOINT` and the bucket name) are constants under `assets/`.
  Every deCODE task holds the `S3ObjectRef`s it reads: the template and per-protein tasks one
  each, the variant info task its file, and `DecodeWindowReader` one per protein.
- **New WF capability** `s3_access: S3Access`, an abstract base class in `build_system/wf/`
  (like `WFDownloader` and `SynapseDownloader`), defaulted in `make_wf`. The production class and
  the test fake both subclass it, so an IDE's go-to-base leads from either one to the interface. It is generic and knows nothing about deCODE:
  - `head(ref: S3ObjectRef) -> S3ObjectInfo` (size, and the md5 from `x-amz-meta-md5chksum`
    when present);
  - `presigned_url(ref: S3ObjectRef, expires_s: int) -> str`.

  The production implementation creates one `boto3` S3 client per
  `(endpoint url, region, profile)` the first time it sees one (SigV4, path-style addressing),
  caches it behind a lock, and reuses it. Tasks get presigned URLs only at execution time.
  Tests inject a fake `S3Access` that maps refs to local files.
- **Rotating keys mid-run.** A cached client keeps the credentials it was created with. On a 403,
  `S3Access` drops that client, rebuilds it once (re-reading the profile), and retries the call.
  So keys renewed with `aws configure` during a long build take effect without a restart.
- **Failure messages.** A missing profile, or a 403 that persists after the rebuild, fails fast.
  The message names the profile and the `aws configure --profile …` command.
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

**Assay set and annotation from Eldjarn et al. 2023, table ST29.** Table ST29
(`ST29_protein_classification` in the supplementary workbook
`41586_2023_6563_MOESM3_ESM.xlsx`) gives, for each assay:
- `SeqId`, in underscore form (`10001_7`);
- `UniProt` (comma-separated for multi-protein targets) and `Gene`;
- `Targeted Soma`, `cis pQTL soma` and `confidence tier`.

It lists exactly the paper's 4,907 analysed SomaScan assays, so it serves both as the **filter**
(the 377 other files are not built) and as the **UniProt source**. A small committed CSV extract
is kept with its provenance in the module docstring: the paper DOI, the workbook file name and
the sheet name. That extract has SeqId, UniProt, Gene, the cis pQTL flags for both versions
where the table gives them, and the confidence tier.

ST2 (`ST2_soma_assays`, 5,032 rows) has fuller assay details. The builder checks it for an
Entrez id column. Otherwise the gene-coordinate reader keys on UniProt and the gene symbol, as
the existing readers allow (§10). The CSF manifest (`csf_aptamer_manifest.csv`, SomaScan 7k,
4,813 of the seq ids) is a cross-check, not a source. The builder records disagreements on
UniProt or gene and does not resolve them silently.

**Generating the manifests.** A one-off regeneration script (like `regenerate_ppp_manifest.py`):
1. lists both S3 prefixes;
2. keeps the files whose seq id is in the ST29 extract;
3. sends a HEAD to each kept object (for size and md5);
4. parses the file names;
5. writes one committed CSV per group under
   `mecfs_bio/assets/reference_data/decode_somascan_sumstats/`.

Columns:
- `seq_id`, in dash form (`10000-28`, from `Proteomics_*_PC0_10000_28_…`), matching the CSF
  convention;
- `gene_symbol` and `target_name` from the file name;
- `key`, `size`, `md5`;
- `uniprot` from ST29.

The script asserts there are 4,907 rows per group, each seq id appearing once.

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

Inputs: the group's index task and the manifest row, which gives the protein's `S3ObjectRef`,
size and md5. In `execute`:
1. `wf.s3_access.head(ref)` and check the size and md5 against the manifest. Then get
   `wf.s3_access.presigned_url(ref, …)`, download with parallel ranges into `scratch_dir`, and
   check the md5.
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
own task is rejected because it duplicates ~39 GB.

### 7. Asset generators and asset layout

- **Generator:** `generate_decode_somascan_protein_tasks(group, index_task, index_name)` returns
  a collection of per-protein tasks, one per manifest row. `group` is a
  `Literal["smp", "raw"]`.
- **Asset ids and paths** are namespaced by group and index, e.g.
  `decode_somascan_smp_hapmap_3_<GENE>_<seq_id_token>`.
- **Modules:** `mecfs_bio/assets/gwas/decode_somascan/` gets `raw/` (manifest-driven
  specifics), `decode_somascan_database/hapmap3/` (index, protein files, hit table), and a coloc
  asset-generator group mirroring CSF's.
- Both groups are defined. The SMP group is built first. The non-normalized group's tasks
  exist but are not built until the sensitivity analysis is wanted.

### 8. Hit table

`ProteinHitTableTask.create(...)` over each group's per-protein tasks, with
`protein_id = ProteinId(seq_id)` and `gene = gene_symbol`, at p < 5e-8, like the CSF hit table
(`hapmap3_csf_hit_table.py`).

### 9. `DecodeWindowReader` (coloc)

Implements `ProteinWindowReader`. It is built with:
- the group's per-protein tasks (to find each protein's seek index and position map in the local
  asset store, via `fetch`);
- the `S3ObjectRef` per protein (from the manifest);
- the variant info asset (for EAF).

Credentials are reached only through `wf.s3_access` at fetch time. The reader holds no keys or
URLs.

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
- **An asset-generator group** (like CSF's) producing the three coloc outputs for a trait (raw
  results, PP.H4 display table, stacked locus plots), plus a `mecfs_bio/analysis/` driver module
  with one `run_*` function per trait. This design adds one run, the exploratory DecodeME run
  (§12).

### 11. LDSC heritability (SMP group)

All-variants and cis-excluded LDSC h² for every protein, as `PppProteinHeritabilityTask`
produces for UKB-PPP (`eur_discovery_hapmap3_ppp_heritability.py`).

- **Reuse the computation, generalize the inputs.** The batched kernel (`batched_ldsc_h2.py`),
  the shared regression context (`PppLdscContext`) and the cis-exclusion logic are
  source-independent. The PPP task's inputs are not:
  - its proteins are `BuildSlimProteinParquetTask`s keyed by Olink id;
  - its N comes from the PPP sample-size table or a constant-N slim column;
  - its gene coordinates come from Sun et al. ST3.

  The builder generalizes it the way `HitTableProtein` already did for the hit table: a source-
  neutral protein entry with the slim task, `ProteinId` and gene, plus a gene-coordinate mapping
  `ProteinId -> GenomicInterval` passed in. The PPP asset (ids, paths, numbers) must come out
  unchanged. If generalizing in place is too invasive, a sibling deCODE task that calls the same
  kernel and context functions is acceptable.
- **Sample size (my call).** The kernel takes one N per protein (`n: (K,)`), but deCODE's N
  varies per variant (about 30.4k–35.7k, roughly ±8%). Use each protein's **median N over its
  regression SNPs**, so the validated kernel stays unchanged. The builder reports the per-protein
  spread of N, for example the 5th to 95th percentile ratio. If it is much wider than in the
  profiled file, they stop and raise it, rather than extending the kernel to per-SNP N.
- **Gene coordinates for cis exclusion** come from the deCODE gene-coordinate reader (§10).
  Proteins without coordinates get all-variants h² only, as PPP handles missing coordinates.
- **LD scores and population:** the same reference as PPP
  (`THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE`, 1000G EUR). Icelanders are European,
  but they are not the reference's population, which is a modest, standard approximation.
- **Output:** the same table shape as the PPP heritability asset, under a deCODE asset id
  (e.g. `decode_somascan_smp_heritability_hapmap_3`).

### 12. Exploratory DecodeME × deCODE coloc run

One `TraitProteinColocTask` run of DecodeME GWAS-1 against the SMP group, in the style of
`decode_me_csf_coloc_exploratory.py`:
- `ColocRunConfig.lead_p_threshold = 1e-5`. `exclude_mhc` follows the existing exploratory
  DecodeME runs (`False` in the CSF exploratory run).
- **Budgets never bind (user decision).** Windows are cheap with the seek index (~1.5 s and
  ~6 MB each), so the run tests every trait-locus × protein candidate pair the hit table
  produces. `SelectionConfig.max_pairs` and `max_per_locus` are set well above the number of
  candidates, for example 10⁶; both must be ≥ 1 integers. The builder adds a check, a log line
  plus an assertion in the run or a column in the outputs, that **no candidate was dropped for
  budget**. A number that only looks big enough is not proof.
- **Expected size.** SomaScan has trans pQTLs for 96% of assays (Table 1), and a DecodeME locus
  near a highly pleiotropic variant can pull in hundreds of proteins. The builder logs candidate
  counts per locus. Run time stays reasonable even at a few thousand pairs: proteins are fetched
  in parallel, and each needs one HEAD plus ~6 MB per window.
- **Interpretation:** pairs at ST36 pleiotropic variants are flagged (see Risks). The module
  docstring labels the run exploratory.
- **Asset ids** are namespaced as the existing exploratory runs are (e.g.
  `decode_me_gwas_1_lead_p_1e-05_x_decode_somascan_smp`, following the established pattern),
  with a `run_*` function in the driver.

## Performance and storage budget (per SomaScan group)

| | Per protein | 4,907 assays |
|---|---|---|
| Build time (download-bound) | ~54 s | ~3.1 days |
| Slim parquet | 7.9 MB | ~39 GB |
| Seek index (zstd) | 0.92 MB | ~4.5 GB |
| Position map | 0.03 MB | ~0.15 GB |
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
- **`S3Access`:** client reuse per (endpoint, region, profile), and the rebuild-once-on-403 path,
  both with an injected fake client factory.
- **Directory-asset dataframe read** (§6): a small build-system test.
- **Heritability (§11):** if the PPP task is generalized, the existing PPP heritability tests must
  pass unchanged. A deCODE-shaped synthetic test covers per-variant N and the median-N rule.
- **Budget check (§12):** a unit test in which a candidate would be dropped makes the
  no-drop check fail.
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
- **SMP artifacts** (see "Choice of SomaScan normalization").
  - Trans pQTLs at highly pleiotropic variants can be reshaped or sign-flipped by the
    normalization. The coloc results already record cis status. The builder also flags coloc
    pairs whose protein lead sits at a variant listed in ST36 (highly pleiotropic in SMP or PC0),
    for example an `is_pleiotropic_variant` column. Mitigation is follow-up analysis, not
    filtering.
  - Genetic correlation with traits tied to overall protein levels (BMI, kidney, liver,
    inflammation) is distorted in both versions, in opposite ways. Interpret it with the
    non-normalized group as a sensitivity check.
- **Comparing plasma and CSF results.** How Western et al. normalized their SomaScan 7k data
  hasn't been checked. Check before comparing deCODE (SMP) and CSF effect sizes or
  correlations directly.

## Left to the builder

- Exact defaults: download chunk size and connection count, `rapidgzip` workers and chunk size,
  `DECODE_MAX_PARALLEL_FETCHES`.
- Reducing memory, against the target above.
- The exact form of the directory-member dataframe read (§6), kept minimal and generic.
- Module and class names beyond those given here.

## Decisions

- **2026-11-21 (user):** SomaScan only; build the SMP group first, with non-normalized as a later
  sensitivity build; database restricted to the paper's 4,907 assays, annotated from ST29.
- **2026-11-21 (user):** S3 access follows option A. Tasks hold an `S3ObjectRef` (endpoint with
  credentials profile name, bucket, key); `WF` provides a generic `S3Access` with no
  bucket or endpoint of its own.
- **2026-11-21 (user):** scope includes LDSC heritability (§11) and an exploratory DecodeME coloc
  run at lead p < 1e-5 with budgets that never bind (§12). rg, LCV and other traits are
  follow-ups.

## Calls made by Claude (for the user to check)

- N for LDSC is each protein's median N over its regression SNPs (§11).
- The profile name `decode-download` is a task-level constant taken from deCODE's instructions
  (§1).
- Directory-asset output with a generic member-file read for slim consumers (§6).
- EAF from `variant_info_iceland`, gated on an all-rows `A1 == effectAllele` check (§4).
- Peak-memory target ≤ 2 GB per protein task (§6).
- Pleiotropic-variant flag (ST36) on coloc results (Risks).
