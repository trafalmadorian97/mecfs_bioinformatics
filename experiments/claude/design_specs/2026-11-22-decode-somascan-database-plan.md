# deCODE SomaScan Database Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A per-group (SMP, non-normalized) deCODE Iceland SomaScan proteomics database with a HapMap3
index, slim per-protein files, gzip seek indexes and position maps. It feeds the hit table,
LDSC h², and a `DecodeWindowReader` that reads remote coloc windows. One exploratory
DecodeME × deCODE coloc run uses it.

**Architecture:** Tasks hold `S3ObjectRef`s (endpoint, bucket, key, AWS profile name). `WF`
provides generic `S3Access` and a parallel range downloader. There is one serial,
resumable task per protein. Each task downloads the 970 MB file with ~128 range requests and
makes one `rapidgzip` pass that yields the HapMap3 rows, a position map and a GZIDX seek
index. Then it deletes the download. Coloc reads windows remotely through `indexed_gzip` and
range requests. EAF comes from `variant_info_iceland.txt.gz` through one rule: true EAF or null.

**Tech Stack:** Python 3.13, polars, numpy, attrs, boto3, rapidgzip 0.16.2, indexed_gzip 1.10.3,
zstandard 0.25.0, requests, pytest.

**Spec:** `experiments/claude/design_specs/2026-11-21-decode-somascan-database-design.md`. Read
it first; § references below point into it. Working code for most algorithms exists in
`experiments/claude/decode_proteomics_survey/mockup_single_protein.py` (functions named below).
Port it into production style; don't re-derive it.

## Global Constraints

- Branch `decode-somascan-database`. **Before Task 1, rebase it onto `origin/main`**, which now
  contains PR #1416 (`ProteinWindowReader` and `DenseWindowFetcher` are ABCs).
- New dependencies go in the pixi **analysis feature** (`feature.analysis.dependencies`,
  conda-forge): `boto3`, `rapidgzip ==0.16.2`, `indexed_gzip ==1.10.3`, `zstandard ==0.25.0`.
- Credentials: AWS profile `decode-download`, region `us-east-1`, endpoint
  `https://s3-ext.decode.is`, bucket `largescaleplasma-2023`. Never commit keys or presigned
  URLs, and never put either in a Task field. Constants live under `mecfs_bio/assets/`; classes live in
  `mecfs_bio/build_system/`.
- Conventions:
  - attrs `@frozen`;
  - ABCs over Protocols for new interfaces;
  - `Literal` / `NewType` ids;
  - `call_with_retries` is the only backoff loop;
  - polars;
  - imports at top of file;
  - no monkeypatching, so inject fakes;
  - Task-level tests preferred;
  - docstrings without backticks or RST;
  - column-name constants, not repeated literals.
- **Columns mean what they say.** EAF is the effect allele's frequency or null, never ImpMAF.
- Group literal: `DecodeSomascanGroup = Literal["smp", "raw"]`. S3 prefixes are
  `final_somascan_smp/` and `final_somascan_raw/`.
- Chromosome code: `chr1`..`chr22` → 1..22, `chrX` → 23. deCODE files list chromosomes in
  **text order** (chr1, chr10, …, chr9, chrX). Window searches go by file order.
- deCODE columns (exact header): `Chrom Pos Name rsids effectAllele otherAllele Beta Pval
  minus_log10_pval SE N ImpMAF`.
- Speed first. Memory only has to stay clear of OOM on the 16 GB machine (4.2 GB per protein is
  fine). Never trade speed for memory (no fewer rapidgzip workers, no smaller blocks).
- After each task: `pixi r invoke green 2>&1 | tee <scratchpad>/green.log`. Check the exit code and
  the log, because the tail alone proves nothing. Stage explicit paths only, and never commit logs.
- Unit tests are offline and fast, about 1 s or less each. Size synthetic data by the seek
  points and chunks a test needs, not by realism: rapidgzip `chunk_size` 32 KiB gives ~30
  seek points per 2 MB. Build expensive fixtures once per module. In targeted pytest runs, list explicit unit-test paths and never the
  `test_mecfs_bio` root.

## Review Focus

1. **Shuffled rows or a different row count in one protein file.** The per-protein task must fail
   loudly, not write misaligned slim rows (Task 9 tests a swapped row and a wrong total).
2. **Window at a chromosome boundary in text order**, e.g. the last rows of chr19 followed by
   chr2. The reader must not read chr2 rows into a chr19 window, or stop early (Task 12:
   `test_window_at_text_order_chromosome_boundary`).
3. **Keys rotated or expired mid-run.** One 403 rebuilds the client and retries. A second 403 fails
   with a message naming `aws configure --profile decode-download` (Task 2).
4. **Null EAF in a window reaching orientation.** A palindrome with null EAF is dropped and counted.
   A non-palindromic SNV is oriented by FASTA (Task 12:
   `test_null_eaf_rows_survive_orientation`).
5. **Remote object changed since the build** (same key, new bytes). The window reader fails with a
   non-retryable error, not wrong data (Task 12: `test_md5_mismatch_is_not_retried`).

---

### Task 1: Dependencies and gzip seek-index compatibility

**Files:**
- Modify: `pyproject.toml` (pixi analysis feature), `pixi.lock`
- Create: `mecfs_bio/util/file_io/gzip_seek_index.py`
- Test: `test_mecfs_bio/unit/util/file_io/test_gzip_seek_index.py`

**Interfaces:**
- Produces:
  - `export_seek_index(gz_path: Path, out_path: Path, parallelization: int, chunk_size: int) -> None`
    makes a rapidgzip pass that writes a zstd-compressed GZIDX v1 (`IndexFormat.INDEXED_GZIP`). Use
    it only in tests. Task 9 exports from its own pass via `export_index_from_open_file` below.
  - `export_index_from_open_file(f: rapidgzip.RapidgzipFile) -> bytes` returns the zstd-compressed GZIDX.
  - `load_seek_index(blob: bytes) -> io.BytesIO` decompresses zstd into a file object for
    `indexed_gzip.import_index`.
  - `SeekPoint(uncompressed: int, compressed: int)` (frozen).
  - `seek_points(blob: bytes, source: BinaryIO) -> tuple[SeekPoint, ...]` imports into an
    `IndexedGzipFile(fileobj=source, auto_build=False)` and reads `seek_points()`.

- [ ] **Step 1:** Add the four dependencies with `pixi add --feature analysis …` at the pinned
  versions, then run `pixi r python -c "import boto3, rapidgzip, indexed_gzip, zstandard"`.
  Expected: no output.
- [ ] **Step 2: Write the failing test** `test_rapidgzip_export_imports_into_indexed_gzip`:
  - write a ~2 MB synthetic gzip of random tab-separated text lines with plain `gzip`
    (`compresslevel=1`);
  - export with `parallelization=4, chunk_size=32 KiB`. Small chunks give ~30 seek points
    from 2 MB, so the test stays small;
  - for 50 random `(offset, length)` pairs, seek and read through `IndexedGzipFile` with the
    imported index;
  - `assert` the bytes equal `data[offset:offset+length]`;
  - `assert len(points) > 10`.

  Measured 2026-11-23: 0.13 s. A 40 MB file with 1 MiB chunks took 5–6 s for the same check.
- [ ] **Step 3:** Run `pixi r pytest test_mecfs_bio/unit/util/file_io/test_gzip_seek_index.py -v`.
  Expected: FAIL (module missing).
- [ ] **Step 4:** Implement the functions, porting `one_pass` (export part) and the import code
  from `mockup_single_protein.py`.
- [ ] **Step 5:** Run the test. Expected: PASS. Then run green.
- [ ] **Step 6:** Commit `feat: rapidgzip/indexed_gzip seek index helpers + compatibility test`.

### Task 2: S3 value types and the `S3Access` WF capability

**Files:**
- Create: `mecfs_bio/build_system/wf/s3_access.py`
- Modify: `mecfs_bio/build_system/wf/base_wf.py`
- Test: `test_mecfs_bio/unit/build_system/wf/test_s3_access.py`

**Interfaces:**
- Produces, all in `s3_access.py`:
  - `S3Endpoint(url: str, region: str, credentials_profile: str)` (frozen).
  - `S3ObjectRef(endpoint: S3Endpoint, bucket: str, key: str)` (frozen).
  - `S3ObjectInfo(size: int, md5: str | None, etag: str)`. `md5` is lowercase hex decoded
    from base64 `x-amz-meta-md5chksum`, or None. `etag` has its quotes stripped.
  - `class S3Access(ABC)` with abstract `head(ref) -> S3ObjectInfo` and
    `presigned_url(ref, expires_s: int) -> str`.
  - `class Boto3S3Access(S3Access)`, constructed as
    `Boto3S3Access(client_factory: S3ClientFactory = Boto3ClientFactory())`.
  - `S3ClientFactory(ABC).make(endpoint: S3Endpoint) -> S3Client`, where `S3Client` is the
    small ABC this module needs: `head_object(bucket, key) -> dict` and
    `generate_presigned_url(bucket, key, expires_s) -> str`.
  - `Boto3ClientFactory` builds `boto3.Session(profile_name=…).client("s3", endpoint_url=…,
    region_name=…, config=Config(signature_version="s3v4", s3={"addressing_style": "path"}))`.
  - `class S3CredentialsError(Exception)`.
  - `WF.s3_access: S3Access`; `make_wf(..., s3_access: S3Access | None = None)` defaults to
    `Boto3S3Access()`.
- Behaviour:
  - Clients are cached per `S3Endpoint` behind a `threading.Lock`.
  - A 403 (`botocore` `ClientError` with HTTP 403) drops that cached client, rebuilds it once
    and retries once.
  - A second 403, or a missing profile (`ProfileNotFound`), raises `S3CredentialsError`. Its
    message comes from the module constant `S3_CREDENTIALS_HELP`, which names the profile and
    `aws configure --profile <profile>`.

- [ ] **Step 1: Write the failing tests** with a fake `S3ClientFactory` that counts `make` calls and
  returns scripted fake clients:
  - `test_client_reused_per_endpoint`: two `head` calls on the same endpoint give `make` count 1.
    A second endpoint gives count 2.
  - `test_403_rebuilds_client_once_then_succeeds`: the first client raises 403 and the second
    returns metadata, so the result arrives and `make` count is 2.
  - `test_persistent_403_raises_credentials_error`: `pytest.raises(S3CredentialsError)`.
  - `test_head_decodes_md5_metadata`: the base64 of bytes `0x93 0x3a…` gives hex `md5`.
    Without the metadata, `md5 is None`, and `etag` has no quotes.
- [ ] **Step 2:** Run `pixi r pytest test_mecfs_bio/unit/build_system/wf/test_s3_access.py -v`.
  Expected: FAIL.
- [ ] **Step 3:** Implement, and wire `make_wf`.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: generic S3Access WF capability with profile credentials`.

### Task 3: Parallel range downloader

**Files:**
- Create: `mecfs_bio/util/download/parallel_ranges.py`, `mecfs_bio/build_system/wf/parallel_range_downloader.py`
- Modify: `mecfs_bio/build_system/wf/base_wf.py`
- Test: `test_mecfs_bio/unit/util/download/test_parallel_ranges.py`

**Interfaces:**
- Produces:
  - `RangeGetter(ABC).get(url: str, start: int, end_inclusive: int) -> RangeResponse`, where
    `RangeResponse(status: int, body: bytes)`. The production version is
    `RequestsRangeGetter` (thread-local `requests.Session`, timeout 120 s).
  - `download_with_parallel_ranges(url: str, size: int, out_path: Path, n_connections: int = 128,
    chunk_bytes: int = 4 * MIB, getter: RangeGetter = RequestsRangeGetter()) -> None`:
    - preallocates with `os.truncate`;
    - runs a thread pool of `n_connections`;
    - `os.pwrite`s each chunk at its offset;
    - wraps each chunk in `call_with_retries`, retrying `requests.ConnectionError`,
      `requests.Timeout`, `ChunkedEncodingError` and a 5xx (raised as a private retryable error);
    - asserts status 206 and an exact body length.
  - `ParallelRangeDownloader(ABC).download(url, size, local_path, md5_hash: str, n_connections,
    chunk_bytes) -> None` downloads and then checks the md5. `HttpParallelRangeDownloader` is the
    production class. `WF.parallel_range_downloader` is wired in `make_wf`, with
    `WF.download_url_parallel(...)` delegating to it.
  - `md5_of_file(path: Path) -> str` lives in the util module; use it, and don't duplicate it in later tasks.

- [ ] **Step 1: Write the failing tests** with a fake `RangeGetter` over an in-memory 10 MB byte string:
  - `test_reassembles_exact_bytes`: with `n_connections=8, chunk_bytes=1 MiB+7`, the output
    equals the source.
  - `test_retries_transient_5xx`: the getter returns 503 once per chunk, then 206, and the
    output is still exact.
  - `test_rejects_200_full_body`: the getter returns 200, which fails with `AssertionError`.
  - `test_md5_mismatch_fails`: `HttpParallelRangeDownloader` with the fake getter and a wrong md5
    raises.
- [ ] **Step 2:** Run the tests. Expected: FAIL.
- [ ] **Step 3:** Implement, porting `parallel_download` / `fetch_range` / `chunk_bounds` from the
  mockup.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: parallel HTTP range downloader as a WF capability`.

### Task 4: Dataframe read of a file inside a `DirectoryAsset` (§6)

**Files:**
- Modify:
  - `mecfs_bio/build_system/meta/read_spec/dataframe_read_spec.py`;
  - `mecfs_bio/build_system/meta/base_meta.py`;
  - `mecfs_bio/build_system/meta/read_spec/read_dataframe.py`;
  - `mecfs_bio/build_system/task/task_util.py`;
  - `mecfs_bio/build_system/rebuilder/metadata_to_path/simple_meta_to_path.py`.
- Create: `mecfs_bio/build_system/meta/gwas_summary_directory_meta.py`
- Test: `test_mecfs_bio/unit/build_system/meta/test_directory_member_read.py`

**Interfaces:**
- Produces:
  - `DataFrameReadSpec.directory_member: PurePath | None = None`. It names the file inside a
    directory asset, and is None for file assets.
  - `DirMeta.read_spec` is a non-abstract property that returns None, like `FileMeta`.
  - `GWASSummaryDataDirectoryMeta(DirMeta)` is frozen, with fields
    `id: AssetId, trait: str, project: str, sub_dir: str, project_path: PurePath,
    read_spec: DataFrameReadSpec | None`. Its path is resolved like `GWASSummaryDataFileMeta`'s,
    minus the extension. Add the dispatch branch next to that one.
  - `scan_dataframe_asset` also accepts `DirectoryAsset` + `DirMeta`. It asserts
    `read_spec.directory_member is not None` and reads `asset.path / directory_member`.
    A `FileAsset` with `directory_member` set is an assertion error.
  - `produces_dataframe` is also true for a `DirMeta` whose `read_spec` is a `DataFrameReadSpec`.

- [ ] **Step 1: Write the failing tests:**
  - `test_scan_reads_member_parquet_of_directory_asset`: write `d/slim.parquet` and scan it
    via `DirectoryAsset(d)` and the meta with `directory_member=PurePath("slim.parquet")`. The
    frame equals what was written.
  - `test_produces_dataframe_for_directory_meta_with_member`.
  - `test_file_asset_with_member_is_rejected`: `pytest.raises(AssertionError)`.
- [ ] **Step 2:** Run the tests. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the tests and the existing meta/read-spec tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: dataframe reads of a named member file of a directory asset`.

### Task 5: deCODE format helpers and the EAF rule (§4)

**Files:**
- Create:
  - `mecfs_bio/constants/decode_somascan_constants.py` (`DecodeSomascanGroup`, column
    constants, `DECODE_HEADER`, `DECODE_TOTAL_ROWS_METADATA_KEY = "total_rows"`, and
    `EAF_TOLERANCE = 1e-5`);
  - `mecfs_bio/build_system/task/decode_database/__init__.py`;
  - `mecfs_bio/build_system/task/decode_database/decode_sumstats_format.py`;
  - `mecfs_bio/build_system/task/decode_database/effect_allele_frequency.py`.
- Test: `test_mecfs_bio/unit/build_system/task/decode_database/test_effect_allele_frequency.py`

**Interfaces:**
- Produces:
  - `chrom_code(chrom: str) -> int` and `chrom_code_expr(col: str) -> pl.Expr` (X → 23). Both
    assert on anything else.
  - `parse_decode_lines(text: bytes) -> pl.DataFrame` gives the full deCODE schema with no
    header (ported from the mockup's `parse_lines` and `FULL_SCHEMA`).
  - `check_header(line: bytes) -> None` asserts the header equals `DECODE_HEADER`.
  - `resolve_effect_allele_frequency(rows: pl.DataFrame, variant_info: pl.DataFrame) -> pl.Series`:
    - `rows` has `CHR` (Int8 code), `POS`, `EA` and `ImpMAF`;
    - `variant_info` has the deduplicated `CHR`, `POS`, `A1` and `ImpFreqA1`;
    - it returns Float64 EAF aligned to `rows` order, named `EAF` (`GWASLAB_EFFECT_ALLELE_FREQ_COL`);
    - the rule: candidates are variant_info rows on (CHR, POS, A1 == EA). A candidate agrees when
      `|min(f, 1-f) - ImpMAF| <= EAF_TOLERANCE`. EAF = f when at least one candidate agrees and
      the agreeing f's span `<= EAF_TOLERANCE`, else null.
    - It never returns ImpMAF.

- [ ] **Step 1: Write the failing test** `test_eaf_rule_cases`. One `rows` frame has one row per
  case, and the assertions are:
  - single agreeing candidate (f=0.7, ImpMAF=0.3) gives 0.7;
  - several candidates where one agrees (f ∈ {0.2, 0.05}, ImpMAF 0.05) gives 0.05;
  - several agreeing with equal f (two A1 rows, both 0.4, ImpMAF 0.4) gives 0.4;
  - f and 1−f both agree (0.3 and 0.7, ImpMAF 0.3) gives null;
  - no candidate gives null;
  - a candidate that disagrees (f=0.2, ImpMAF 0.4) gives null.

  Also `test_eaf_preserves_row_order` (shuffled input gives the matching output order) and
  `test_chrom_code_text_order_and_x`.
- [ ] **Step 2:** Run `pixi r pytest test_mecfs_bio/unit/build_system/task/decode_database/test_effect_allele_frequency.py -v`.
  Expected: FAIL.
- [ ] **Step 3:** Implement. Add a row index, left-join, group by row index, then sort back; the
  logic is that of `check_variant_info_maf_disambiguation.py`.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: deCODE format helpers and the true-EAF-or-null rule`.

### Task 6: Synthetic deCODE fixtures and `BuildDecodeVariantInfoTask` (§4)

**Files:**
- Create:
  - `test_mecfs_bio/unit/build_system/task/decode_database/decode_fixtures.py`;
  - `mecfs_bio/build_system/task/decode_database/build_decode_variant_info_task.py`.
- Test: `test_mecfs_bio/unit/build_system/task/decode_database/test_build_decode_variant_info_task.py`

**Interfaces:**
- Produces:
  - Fixtures, used by every later deCODE test:
    - `write_synthetic_decode_file(path, rows: pl.DataFrame) -> None` writes the exact deCODE
      text with header as plain gzip, with `compresslevel=1` so that blocks span seek points.
    - `synthetic_rows(seed: int) -> pl.DataFrame` gives ~20k rows (~1.3 MB of text) over chr1,
      chr10, chr2, chr22 and chrX in text order. It includes one `!` row, one multi-allelic site with f/1−f, and
      positions dense enough for several windows.
    - `write_synthetic_variant_info(path, rows) -> None` writes 5 columns, each row written
      twice (the duplicates), plus multi-allelic rows.
    - `FakeS3Access(S3Access)` maps `S3ObjectRef` to a local path. Its `head` returns the size
      and md5, and `md5=None` with the md5 as etag when built with `etag_only=True`. Its
      `presigned_url` returns `file://` paths.
    - `FakeParallelRangeDownloader` copies `file://` sources and checks the md5.
    - `CountingRangeSource` serves byte ranges of a local file and counts bytes.
  - `BuildDecodeVariantInfoTask(GeneratingTask)` has fields `meta`, `ref: S3ObjectRef` and
    `n_connections: int = 128`, built by `create(asset_id: str, ref: S3ObjectRef)`.
    - `execute`: `head`; assert md5 is None, then use `etag` as the md5; download with parallel
      ranges; stream into parquet with `pl.scan_csv(..., separator="\t")`; `unique()`; cast
      `CHR`(Int8 code), `POS`(Int64), `A1`(String), `ImpFreqA1`(Float64), `ImpInfo`(Float32);
      sort by (CHR, POS); write zstd with `row_group_size=65_536`.
    - Output is a `FileAsset` with a `ReferenceFileMeta` and a parquet read spec.
    - Column names are constants `VI_CHR_COL`, `VI_POS_COL`, `VI_A1_COL`, `VI_FREQ_COL`,
      `VI_INFO_COL` in `decode_somascan_constants.py`.

- [ ] **Step 1: Write the failing tests:**
  - `test_variant_info_dedups_sorts_and_checks_etag`: the output row count equals the distinct
    rows, it is sorted by (CHR, POS), and every row-group size is ≤ 65,536 (via `pyarrow.parquet.ParquetFile`).
  - `test_variant_info_etag_mismatch_fails`.
- [ ] **Step 2:** Run the tests. Expected: FAIL.
- [ ] **Step 3:** Implement the fixtures and the task.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: deCODE variant info task (dedup, sorted, ETag-verified)`.

### Task 7: Manifests (§3)

**Files:**
- Create:
  - `mecfs_bio/assets/reference_data/decode_somascan_sumstats/__init__.py`;
  - `.../decode_somascan_sumstats/st29_soma_assays.csv` (committed ST29 extract);
  - `.../decode_somascan_sumstats/extract_st29.py`;
  - `.../decode_somascan_sumstats/regenerate_decode_manifest.py`;
  - `.../decode_somascan_sumstats/decode_somascan_smp_manifest.csv` and `..._raw_manifest.csv`;
  - `mecfs_bio/assets/gwas/decode_somascan/raw/decode_s3.py` with `DECODE_S3_ENDPOINT`,
    `DECODE_BUCKET`, `DECODE_PREFIX: dict[DecodeSomascanGroup, str]` and
    `DECODE_VARIANT_INFO_REF`;
  - `mecfs_bio/asset_generator/decode_somascan_manifest.py`.
- Test: `test_mecfs_bio/unit/asset_generator/test_decode_somascan_manifest.py`

**Interfaces:**
- Produces:
  - `DecodeSomascanProtein`, frozen, with fields:
    - `seq_id: SeqId`, in dash form;
    - `gene_symbol: str` and `target_name: str`, from the file name;
    - `uniprot: UniProtId`, from ST29 with commas kept;
    - `entrez_gene_id: EntrezGeneId | None`;
    - `key: str`, `size: int`, `md5: str`.

    Property `ref` → `S3ObjectRef(DECODE_S3_ENDPOINT, DECODE_BUCKET, key)`.
  - `load_decode_somascan_proteins(group: DecodeSomascanGroup) -> tuple[DecodeSomascanProtein, ...]`
    follows manifest order, which is sorted by key.
  - Manifest columns: `seq_id, gene_symbol, target_name, uniprot, entrez_gene_id, key, size, md5`.
- Facts settled while writing this plan:
  - ST2 has **no Entrez column**. Its columns are `seqid uniprot target_name target_full_name
    gene_name dilution …`.
  - So `entrez_gene_id` is resolved in the regenerate script, never by gene symbol:
    UniProt accession(s) → HGNC complete set (`HGNC_UNIPROT_GENES`, materialized with
    `DEFAULT_RUNNER`) → `entrez_id`. Several ids are joined with `|`, as in the CSF manifest.
  - ST29 sheet `ST29_protein_classification` has its header on row 3. Its columns are `UniProt,
    OlinkID, SeqId, Gene, …, Targeted Soma, …, cis pQTL soma, …`. Keep the rows with SeqId
    non-null and `Targeted Soma == "Y"`.

- [ ] **Step 1:** Write `extract_st29.py`. It reads the workbook path given as an argument
  (`41586_2023_6563_MOESM3_ESM.xlsx`) and writes `st29_soma_assays.csv` with `seq_id,
  uniprot, gene, cis_pqtl_soma, confidence_tier`. The provenance (DOI 10.1038/s41586-023-06563-x,
  file name and sheet) goes in the module docstring. Run it on the workbook. A copy is at the
  scratchpad `papers_eldjarn/supp_tables.xlsx`, or download it from the paper. Then
  `assert` 4,907 rows with unique `seq_id`.
- [ ] **Step 2: Write the failing test** `test_load_manifest_rows_and_refs` over the committed
  CSVs:
  - 4,907 rows per group;
  - unique `seq_id`;
  - the same `seq_id` set in both groups;
  - every key under its group prefix;
  - md5 is 32 hex characters.
- [ ] **Step 3:** Write `regenerate_decode_manifest.py`, modelled on `regenerate_csf_manifest.py`:
  1. list both prefixes with boto3 under profile `decode-download`;
  2. filter by ST29;
  3. HEAD each kept object for size and md5;
  4. parse file names. The regex `Proteomics_(?:SMP_)?PC0_(\d+)_(\d+)_([^_]+)_(.+)_\d{8}\.txt\.gz`
     gives seq id `\1-\2`, gene and target;
  5. resolve Entrez through HGNC;
  6. write the CSVs;
  7. assert 4,907 rows per group;
  8. print disagreements between ST29 `Gene` and the file-name gene, and against
     `csf_aptamer_manifest.csv` UniProt/gene. Print them only; don't resolve them.
- [ ] **Step 4:** Ask the user to confirm the `decode-download` profile exists (`aws configure
  --profile decode-download`, region `us-east-1`; keys from their scratchpad `env.sh` are fine).
  Then run `pixi r python mecfs_bio/assets/reference_data/decode_somascan_sumstats/regenerate_decode_manifest.py
  | tee <scratchpad>/manifest.log`. Report the disagreement counts to the user.
- [ ] **Step 5:** Write the loader and run the test. Expected: PASS. Run green.
- [ ] **Step 6:** Commit the CSVs, scripts, loader and constants. Leave the log out:
  `feat: deCODE SomaScan manifests from ST29 and the bucket listing`.

### Task 8: `ConstructDecodeVariantIndexTask` (§5)

**Files:**
- Create: `mecfs_bio/build_system/task/decode_database/construct_decode_variant_index_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/decode_database/test_construct_decode_variant_index_task.py`

**Interfaces:**
- Consumes: Tasks 2, 3, 5 and 6.
- Produces:
  - `ConstructDecodeVariantIndexTask(GeneratingTask)` has fields `meta`,
    `template: DecodeSomascanProtein`, `membership_task: Task` and `variant_info_task: Task`.
    `deps` = [membership_task, variant_info_task]. It is built by
    `create(asset_id: str, template: DecodeSomascanProtein, membership_task, variant_info_task)`.
  - Output parquet columns, in order, using the CSF/PPP constants:
    - `CHR` Int8;
    - `POS` Int64;
    - `EA`, `NEA`, `rsID`;
    - `EAF` Float32;
    - `is_strand_ambiguous`, as the PPP LDSC context requires (`PPP_INDEX_IS_STRAND_AMBIGUOUS_COL`);
    - `SOURCE_ROW` UInt64, with constant `DECODE_INDEX_SOURCE_ROW_COL`;
    - `NAME` String, with constant `DECODE_INDEX_NAME_COL`.

    Rows are in file order (ascending `SOURCE_ROW`). The parquet key-value metadata
    `total_rows` holds the template's data row count.
  - `read_index_total_rows(path: Path) -> int`.
- Algorithm:
  1. Download the template the way Task 9 does: `head`, check against the manifest,
     `presigned_url`, then `download_url_parallel`.
  2. Decompress it with `rapidgzip` and parse it fully with polars.
  3. Drop `otherAllele == "!"` rows.
  4. Inner-join to the membership list on (CHR, POS, unordered allele key), reusing
     `unordered_allele_key_narwhals`. Adopt the deCODE orientation.
  5. Fill EAF with `resolve_effect_allele_frequency` against the variant info rows for the
     index positions.
  6. Assert that EAF is non-null on every row, with a message naming the count and the template.

- [ ] **Step 1: Write the failing tests** on the synthetic fixture, with a small membership frame:
  - `test_index_rows_source_row_and_total_rows`: `SOURCE_ROW` points at the same `Name` in the
    synthetic file, `total_rows` equals the data row count, and the `!` row is absent.
  - `test_index_eaf_is_true_effect_allele_frequency`: a row with EA as the major allele gets
    EAF > 0.5.
  - `test_index_fails_when_eaf_undetermined`: put an f/1−f site in membership, and
    `pytest.raises(AssertionError)`.
- [ ] **Step 2:** Run the tests. Expected: FAIL.
- [ ] **Step 3:** Implement, porting `run_template` from the mockup.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: deCODE shared HapMap3 variant index with SOURCE_ROW and true EAF`.

### Task 9: `BuildDecodeSomascanProteinTask` (§6)

**Files:**
- Create:
  - `mecfs_bio/build_system/task/decode_database/build_decode_somascan_protein_task.py`;
  - `mecfs_bio/build_system/task/decode_database/one_pass.py`.
- Test: `test_mecfs_bio/unit/build_system/task/decode_database/test_build_decode_somascan_protein_task.py`

**Interfaces:**
- Consumes: Tasks 1–5 and 8, and `write_byte_stream_split_parquet`.
- Produces:
  - File-name constants `DECODE_SLIM_FILENAME = "slim.parquet"`,
    `DECODE_SEEK_INDEX_FILENAME = "seek_index.gzidx.zst"` and
    `DECODE_POSMAP_FILENAME = "posmap.parquet"`.
  - Posmap columns: `chrom` Int8, `pos` Int64, `line_offset` UInt64, `row_number` UInt64.
    Posmap metadata keys: `source_key`, `compressed_size`, `md5`, `header_end_offset`,
    `posmap_format_version` = "1" and `gzidx_format_version` = "1".
  - `OnePassResult(slim: pl.DataFrame, posmap: pl.DataFrame, header_end_offset: int, seek_index: bytes)`.
  - `one_pass(gz_path: Path, index: pl.DataFrame, total_rows: int, parallelization: int = 8,
    chunk_size: int = 16 * MIB, posmap_step: int = 1 * MIB) -> OnePassResult` asserts:
    - the header;
    - the total row count;
    - gathered `Name` == index `NAME`, row for row;
    - no gathered `!`.

    Port `line_blocks`, `posmap_rows`, `gather_lines` and `one_pass` from the mockup.
  - `BuildDecodeSomascanProteinTask(GeneratingTask)` has fields `meta: GWASSummaryDataDirectoryMeta`,
    `index_task: Task` and `protein: DecodeSomascanProtein`, with `deps` = [index_task].
    `create(index_task, protein, asset_id, index_name, group)` sets
    `read_spec=DataFrameReadSpec(DataFrameParquetFormat(), directory_member=PurePath(DECODE_SLIM_FILENAME))`,
    `trait="decode_somascan_<group>"` and `project=gene_symbol`.
  - `execute`:
    1. `head` and assert the size and md5 equal the manifest;
    2. `presigned_url(ref, 6 * 3600)`;
    3. `download_url_parallel`;
    4. `one_pass`;
    5. write the three files into `scratch_dir/out`, with slim BETA/SE/N as float32 and
       byte-stream-split zstd;
    6. delete the download;
    7. return `DirectoryAsset`.

    Log the peak RSS (`resource.getrusage`) and the wall time per step.

- [ ] **Step 1: Write the failing tests**, building through the fakes with a small
  `chunk_size` (32 KiB, giving ~20 seek points) and `posmap_step` (8 KiB). Build the outputs
  once in a module-scoped fixture that Task 12's tests also use:
  - `test_protein_outputs_match_brute_force`: the slim BETA/SE/N equal a brute-force parse joined on
    `SOURCE_ROW`. Every posmap row's `line_offset` starts a line whose (chrom, pos) match. The
    seek index round-trips through `seek_points`.
  - `test_protein_fails_on_swapped_rows`: swap two index rows in the protein file, then
    `pytest.raises(AssertionError)`.
  - `test_protein_fails_on_wrong_total_rows`: drop one non-index row.
  - `test_protein_fails_on_md5_mismatch`: the manifest md5 differs from `head`.
  - `test_slim_reads_through_scan_dataframe_asset`: Task 4's path on the real output.
- [ ] **Step 2:** Run the tests. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: per-protein deCODE task (slim, seek index, position map)`.

### Task 10: Asset generator, asset modules, hit tables (§7, §8)

**Files:**
- Create:
  - `mecfs_bio/asset_generator/decode_somascan_protein_asset_generator.py`;
  - `mecfs_bio/assets/gwas/decode_somascan/__init__.py` and `raw/__init__.py`;
  - `mecfs_bio/assets/gwas/decode_somascan/raw/decode_variant_info.py`
    (`DECODE_VARIANT_INFO = BuildDecodeVariantInfoTask.create("decode_variant_info_iceland", DECODE_VARIANT_INFO_REF)`);
  - `mecfs_bio/assets/gwas/decode_somascan/decode_somascan_database/hapmap3/` with `__init__.py`,
    `hapmap3_decode_index.py` (`HAPMAP_3_DECODE_SMP_INDEX`, `HAPMAP_3_DECODE_RAW_INDEX`),
    `hapmap3_decode_database_protein_files.py` (`HAPMAP_3_DECODE_SMP_DATABASE`,
    `HAPMAP_3_DECODE_RAW_DATABASE`) and `hapmap3_decode_hit_table.py`
    (`HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8`, `HAPMAP_3_DECODE_RAW_HIT_TABLE_5E8`);
  - `mecfs_bio/analysis/decode_somascan_database_build.py`, modelled on `csf_pqtl_database_build.py`.
    `main()` builds the SMP variant info, index, proteins and hit table.

**Interfaces:**
- Produces:
  - `DecodeSomascanTaskCollection(group, protein_tasks: tuple[BuildDecodeSomascanProteinTask, ...])`,
    with `hit_table_proteins() -> tuple[HitTableProtein, ...]` using `ProteinId(seq_id)` and
    `gene=gene_symbol`.
  - `generate_decode_somascan_protein_tasks(group, index_task, index_name) -> DecodeSomascanTaskCollection`.
  - Asset ids: index `decode_somascan_<group>_hapmap_3_index`, proteins
    `decode_somascan_<group>_hapmap_3_<GENE>_<seq_id with - → _>`, hit table
    `decode_somascan_<group>_hapmap_3_hit_table_5e8`.
  - The template is `load_decode_somascan_proteins(group)[0]`.

- [ ] **Step 1:** Write the modules. There are no structural tests for generators: the import-time
  construction is the check.
- [ ] **Step 2:** Run `pixi r python -c "import mecfs_bio.assets.gwas.decode_somascan.decode_somascan_database.hapmap3.hapmap3_decode_hit_table as m; print(len(m.HAPMAP_3_DECODE_SMP_HIT_TABLE_5E8.proteins))"`.
  Expected: `4907`. Run green.
- [ ] **Step 3:** Commit `feat: deCODE SomaScan asset generator, both groups, hit tables`.

### Task 11: Window fetchers declare dependencies and get `fetch`

The spec says `DecodeWindowReader` finds each protein's seek index and position map "in the
local asset store, via fetch". Today readers receive only `scratch_dir` and `wf`, and
`TraitProteinColocTask.deps` doesn't include the fetcher's tasks. This task adds both, as a
mechanical interface change.

**Files:**
- Modify:
  - `mecfs_bio/build_system/task/protein_coloc/dense_window_fetch.py`;
  - `.../ppp_window_reader.py` and `.../csf_window_reader.py`;
  - `.../trait_protein_coloc_task.py`;
  - their tests in `test_mecfs_bio/unit/build_system/task/protein_coloc/`.

**Interfaces:**
- Produces:
  - `ProteinWindowReader.deps -> list[Task]`, a concrete property returning `[]`.
  - `ProteinWindowReader.fetch_protein(protein_id, windows, scratch_dir, fetch: Fetch, wf)`.
  - `DenseWindowFetcher.deps -> list[Task]`, concrete and returning `[]`.
  - `DenseWindowFetcher.fetch_all(requests, scratch_dir, fetch: Fetch, wf)`.
  - `ParallelWindowFetcher.deps` returns `self.reader.deps`.
  - `TraitProteinColocTask.deps` appends `self.fetcher.deps`, and `execute` passes `fetch` through.

- [ ] **Step 1: Write the failing test** `test_coloc_task_deps_include_fetcher_deps` in
  `test_trait_protein_coloc_task.py`. A fake reader whose `deps` returns a sentinel task leads to
  that sentinel being in `coloc_task.deps`.
- [ ] **Step 2:** Run it. Expected: FAIL.
- [ ] **Step 3:** Implement and update every call site and fake reader.
- [ ] **Step 4:** Run `pixi r pytest test_mecfs_bio/unit/build_system/task/protein_coloc/ -v`.
  Expected: all PASS. Run green.
- [ ] **Step 5:** Commit `refactor: window fetchers declare task deps and receive fetch`.

### Task 12: `DecodeWindowReader` (§9)

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/decode_window_reader.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_decode_window_reader.py`

**Interfaces:**
- Consumes: Tasks 1, 2, 5, 6, 9 and 11.
- Produces:
  - `DecodeWindowReader(ProteinWindowReader)`, frozen, with fields:
    - `proteins: tuple[DecodeProteinEntry, ...]`, where
      `DecodeProteinEntry(task: BuildDecodeSomascanProteinTask)` is looked up by `seq_id`;
    - `variant_info_task: Task`;
    - `range_source_factory: RangeSourceFactory = PresignedRangeSourceFactory()`;
    - `prefetch_connections: int = 16`, `prefetch_piece_bytes: int = 512 * 1024`.

    `deps` = the protein tasks + `variant_info_task`.
  - `RangeSourceFactory(ABC).open(url: str, size: int) -> RangeSource`, where
    `RangeSource(ABC).read_range(start, length) -> bytes`. The production version uses
    `RequestsRangeGetter` from Task 3.
  - `DECODE_MAX_PARALLEL_FETCHES = 6`, a starting point that Task 17 measures and sets.
  - `class RemoteObjectChanged(Exception)` is not an `OSError`, so `ParallelWindowFetcher`
    doesn't retry it.
  - Helpers ported from the mockup: `uncompressed_span(posmap, window, eof)` (text-order aware),
    `compressed_span(points, span, size)` and `RemoteGzipSource(io.RawIOBase)`, which serves
    prefetched bytes and range-reads the rest.
- `fetch_protein`:
  1. load the posmap and seek index via `fetch(task.asset_id)`;
  2. HEAD the object. A size or md5 different from the posmap metadata raises `RemoteObjectChanged`;
  3. for each window, in file order: find the spans, prefetch, `IndexedGzipFile(fileobj=…,
     index_file=…, auto_build=False)`, seek and read, parse, and filter to the window;
  4. convert to `WINDOW_SCHEMA` (CHR code, EA/NEA, BETA, SE, N). EAF comes from
     `resolve_effect_allele_frequency`, with the variant-info rows read by a lazy
     `(CHR == c) & POS.is_between(lo, hi)` filter on the sorted parquet. `!` rows are dropped.

- [ ] **Step 1: Write the failing tests**, using the Task 9 outputs built in a fixture through the
  fakes, a `CountingRangeSource` factory and a brute-force filter of the synthetic rows:
  - `test_windows_match_brute_force`: windows at file start, file end, spanning a seek point, on
    chr22 (short) and on chrX. Each frame equals brute force in `WINDOW_SCHEMA`, with EAF from
    the rule.
  - `test_window_at_text_order_chromosome_boundary`: a window to the end of chr10 contains no chr2
    rows, and a window at the start of chr2 is complete.
  - `test_prefetch_covers_reads`: the counting source records zero non-prefetched bytes.
  - `test_eaf_null_where_rule_undetermined`: the f/1−f site's row has null EAF, and no EAF equals
    ImpMAF where EA is the major allele.
  - `test_md5_mismatch_is_not_retried`: the fake head returns another md5, which raises
    `RemoteObjectChanged`, and a reader wrapped in `ParallelWindowFetcher` raises once and
    doesn't retry.
  - `test_null_eaf_rows_survive_orientation` (in `test_orient_and_align.py`): a window with one
    null-EAF palindrome and one null-EAF non-palindromic SNV goes through the existing orient
    function. The palindrome is dropped and counted, and the SNV is kept and oriented by FASTA. If
    this fails, stop and report: it is a spec assumption (§4).
- [ ] **Step 2:** Run the tests. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run the tests. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: DecodeWindowReader (remote seek-index window reads, true EAF)`.

### Task 13: Coloc budget-never-binds check (§12)

**Files:**
- Modify: `mecfs_bio/build_system/task/protein_coloc/trait_protein_coloc_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_trait_protein_coloc_task.py`

**Interfaces:**
- Produces:
  - `ColocRunConfig.require_no_budget_drop: bool = False`.
  - When true, `execute` asserts after `select_candidates` that no row has `STATUS ==
    "not_selected_budget"`. The assertion message gives the dropped count and the budgets.
  - `execute` always logs the candidate counts per locus (`n_candidates_by_locus`).

- [ ] **Step 1: Write the failing test**
  `test_require_no_budget_drop_fails_when_budget_binds`. Use the existing task test fixture with
  `SelectionConfig(max_pairs=1, max_per_locus=1)` and more than one candidate, then
  `pytest.raises(AssertionError)`. The same fixture with large budgets passes.
- [ ] **Step 2:** Run it. Expected: FAIL.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Run it. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: coloc run option asserting no candidate is dropped for budget`.

### Task 14: deCODE coloc source, gene coordinates and the exploratory DecodeME run (§10, §12)

**Files:**
- Modify:
  - `mecfs_bio/build_system/task/protein_coloc/csf_gene_coords_reader.py`, extracting the
    per-protein rule into a source-neutral function;
  - `mecfs_bio/asset_generator/protein_coloc_asset_generator.py`, extending
    `ProteomicsSourceInfix` with `"decode_somascan_smp"` and `"decode_somascan_raw"`.
- Create:
  - `mecfs_bio/build_system/task/protein_coloc/decode_gene_coords_reader.py`;
  - `mecfs_bio/asset_generator/decode_somascan_coloc_asset_generator.py`;
  - `mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/decode_me_decode_somascan_coloc_exploratory.py`;
  - `mecfs_bio/analysis/decode_somascan_coloc_analysis.py`.
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_decode_gene_coords_reader.py`

**Interfaces:**
- Produces:
  - `GeneKeys(entrez_ids: frozenset[str], symbol: str)` and
    `gene_interval(keys: GeneKeys, genes: MagmaGenes, conflict_distance_bp: int) -> GenomicInterval | None`.
    This is the CSF rule unchanged. `_MagmaGenes` becomes public as `MagmaGenes`, and
    `magma_genes(table)` is exposed. `CsfGeneCoordsReader` calls them, and its tests pass unchanged.
  - `DecodeGeneCoordsReader(proteins: tuple[DecodeSomascanProtein, ...], conflict_distance_bp=…)`
    is keyed by `ProteinId(seq_id)`, with `entrez_gene_id` split on `|`.
  - `DecodeProteinAccessionsReader(proteins: tuple[DecodeSomascanProtein, ...])`, in
    `protein_accessions_reader.py`. Like `CsfGeneCoordsReader`, it takes its data from the
    manifest rows it holds and ignores the table passed in. It returns one row per
    (seq id, UniProt), with ST29's comma-separated UniProt values split.
    `generate_coloc_protein_info_task` stops asserting that the run is PPP-only. It picks
    `St3ProteinAccessionsReader()` for an `St3GeneCoordsReader` run and
    `DecodeProteinAccessionsReader(reader.proteins)` for a `DecodeGeneCoordsReader` run, and it
    asserts on anything else. This is the existing targeted generator; it doesn't run per coloc run.
  - `generate_decode_somascan_coloc_tasks(base_name, trait_task, trait_model, config, group="smp", h4_threshold=0.8) -> ColocTaskGroup`
    builds `ProteomicsColocSource(asset_infix=f"decode_somascan_{group}", hit table,
    MAGMA GRCh38 table, DecodeGeneCoordsReader, ParallelWindowFetcher(DecodeWindowReader(...),
    DECODE_MAX_PARALLEL_FETCHES))` and calls `generate_protein_coloc_tasks`. No plot tasks.
  - `DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5` uses:
    - `base_name="decode_me_gwas_1_lead_p_1e_5"`;
    - `ColocRunConfig(lead_p_threshold=1e-5, exclude_mhc=False,
      selection=SelectionConfig(max_pairs=1_000_000, max_per_locus=1_000_000),
      require_no_budget_drop=True)`;
    - the case-control trait model as in `decode_me_csf_coloc_exploratory.py`.

    Its docstring labels it exploratory.
  - `decode_somascan_coloc_analysis.py` has `run_decode_me_exploratory()` and a `main()` that calls it.

- [ ] **Step 1: Write the failing tests:**
  - `test_decode_gene_coords_uses_entrez_then_symbol_rules`: one protein each for an id hit, a
    far symbol conflict (None), several ids (None), and the symbol fallback.
  - `test_decode_accessions_split_multi_uniprot` (in `test_protein_accessions_reader.py`): a
    protein with `"P1, P2"` gives two rows in `PROTEIN_ACCESSIONS_SCHEMA`.
  - Run the existing `test_csf_gene_coords_reader.py` unchanged.
- [ ] **Step 2:** Run them. Expected: the new test FAILS and the CSF tests PASS.
- [ ] **Step 3:** Implement the extraction, the reader, the generator, the run module and the driver.
- [ ] **Step 4:** Run the tests, plus an import check of the exploratory module. Expected: PASS. Run green.
- [ ] **Step 5:** Commit `feat: deCODE coloc source and exploratory DecodeME run at lead p 1e-5`.

### Task 15: LDSC heritability for the SMP group (§11)

Approach: a **sibling task** that reuses the kernel, context and cis-mask functions, so PPP stays
byte-identical. Generalizing `PppProteinHeritabilityTask` in place would mean re-keying its Oid,
SynID and sample-size paths, which the spec calls acceptable to avoid.

**Files:**
- Create:
  - `mecfs_bio/build_system/task/ppp_ldsc/decode_protein_heritability_task.py`;
  - `mecfs_bio/assets/gwas/decode_somascan/decode_somascan_database/hapmap3/hapmap3_decode_smp_heritability.py`.
- Modify: `mecfs_bio/build_system/task/ppp_ldsc/ppp_protein_heritability_task.py`. Only move
  `_build_context`, `_build_cis_exclude` and `_read_table` to public names in
  `ppp_ldsc/heritability_common.py`, with no behaviour change.
- Test: `test_mecfs_bio/unit/build_system/task/ppp_ldsc/test_decode_protein_heritability_task.py`

**Interfaces:**
- Produces:
  - `DecodeProteinHeritabilityTask(GeneratingTask)` has fields `meta`,
    `proteins: tuple[DecodeSomascanProtein, ...]`, `protein_tasks`, `index_task`,
    `ld_scores_task`, `gene_coords_task` (the MAGMA GRCh38 table) and `config: PppHeritabilityConfig`.
  - Gene coordinates come from `DecodeGeneCoordsReader`.
  - N is the **median of the slim N over the context SNPs** per protein.
  - Output columns are those of the PPP table, except that the id column is `protein_id`
    (constant `DECODE_H2_PROTEIN_ID_COL`) holding the seq id, because these are not Olink ids. It
    adds `n_p05` and `n_p95` columns, which are the 5th and 95th percentiles of N.
  - Asset `HAPMAP_3_DECODE_SMP_HERITABILITY`, id `decode_somascan_smp_heritability_hapmap_3`.

- [ ] **Step 1: Write the failing test**
  `test_decode_h2_uses_median_n_and_matches_ppp_kernel`. Build a synthetic slim with per-variant
  N varying ±8% around 33,000 on the existing CT-LDSC/PPP synthetic recipe. Assert:
  - `n_bar` equals the median;
  - h² equals `batched_h2` called directly with that median;
  - a protein without coordinates has only the all-variants row.
- [ ] **Step 2:** Run it. Expected: FAIL.
- [ ] **Step 3:** Implement. Then run the existing PPP heritability tests unchanged, and the new test.
  Expected: all PASS. Run green.
- [ ] **Step 4:** Commit `feat: LDSC heritability for deCODE SomaScan proteins (median per-protein N)`.

### Task 16: ST36 pleiotropic-variant flag on coloc results (Risks)

ST36 lists Olink lead positions, with the r² and distance to the SMP and PC0 counterparts, but
not the SomaScan variant itself. A flag named "is the pleiotropic variant" would claim more than
the data say. The column is therefore **`NEAR_ST36_PLEIOTROPIC_SMP_250KB`**: true when the
protein lead lies within 250 kb of an ST36 row marked `Highly Pleiotropic Soma SMP == "Y"`,
measuring from the Olink position plus or minus the SMP distance.

**Files:**
- Create:
  - `.../decode_somascan_sumstats/st36_pleiotropic_pqtls.csv`, extracted by adding an
    `--st36` mode to `extract_st29.py`;
  - `mecfs_bio/build_system/task/protein_coloc/coloc_pleiotropy_flag_task.py`;
  - in the exploratory module: `DECODE_ME_DECODE_SOMASCAN_SMP_COLOC_LEAD_P_1E_5_PLEIOTROPY`.
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_pleiotropy_flag_task.py`

**Interfaces:**
- Produces: `ColocPleiotropyFlagTask.create(asset_id, coloc_task, pleiotropic_csv: Path, distance_bp=250_000)`,
  which writes the raw results plus the flag column as a parquet `ResultTableMeta`.
- [ ] **Step 1: Write the failing test** `test_flag_marks_leads_near_st36_smp_rows`. Use
  `coloc_result_fixtures.py` results with leads at 0 kb, 249 kb and 251 kb from an ST36 SMP row,
  and one on another chromosome. The expected flags are `[True, True, False, False]`, and a
  null protein lead gives null.
- [ ] **Step 2:** Run it. Expected: FAIL. Step 3: implement. Step 4: run it, expect PASS, and run green.
- [ ] **Step 5:** Commit `feat: flag coloc pairs whose protein lead is near an ST36 pleiotropic pQTL`.

### Task 17: Real-bucket smoke run and defaults (manual, needs the profile)

- [ ] **Step 1:** Run `DEFAULT_RUNNER` on `DECODE_VARIANT_INFO`, then `HAPMAP_3_DECODE_SMP_INDEX`, then
  the per-protein tasks for AP4M1 (`10076-1`) and the template. Record the wall time and peak
  RSS per protein in `experiments/claude/decode_proteomics_survey/smoke_run.md`, which is
  committed (the log is not). Expected:
  - index 1,197,450 rows;
  - ~54 s per protein;
  - peak below 5 GB;
  - the slim matches the mockup's slim for AP4M1 (`assert_frame_equal`).
- [ ] **Step 2:** Measure download `chunk_bytes` at 2, 4 and 8 MiB on one protein, and keep the
  fastest. Measure `DECODE_MAX_PARALLEL_FETCHES` ∈ {4, 6, 8} on 24 proteins × 1 window. Set the
  constants and record the numbers.
- [ ] **Step 3:** Read 3 windows for AP4M1 through `DecodeWindowReader` and compare them with a
  brute-force scan of the downloaded file (the mockup's `brute_force_window`). Expected: equal.
- [ ] **Step 4:** Report the per-protein N spread (p05/p95 ratio) on AP4M1. If it is much wider than
  30.4k–35.7k, stop and raise it (§11).
- [ ] **Step 5:** Commit the constants and `smoke_run.md`:
  `chore: deCODE smoke-run measurements and tuned defaults`.
- [ ] **Step 6:** Hand over to the user. They start the full SMP build (~3.1 days,
  `pixi r python mecfs_bio/analysis/decode_somascan_database_build.py`), then heritability and
  `run_decode_me_exploratory()`. Neither is run by the implementer.
