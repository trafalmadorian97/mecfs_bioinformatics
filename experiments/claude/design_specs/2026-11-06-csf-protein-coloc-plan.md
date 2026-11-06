# CSF Protein Colocalization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run coloc.abf between eight trait GWASs (plus two exploratory DecodeME runs) and the 7,008 Western et al. 2024 CSF aptamer GWASs, producing the same raw-results / display-table / locus-plot outputs as the UKB-PPP coloc.

**Architecture:** Bring the CSF database (slim HapMap3 per-aptamer files, already built) from the `csf-database` branch; generalize the PPP-specific pieces of `protein_coloc` (protein id type, gene-coordinate reader, hit table, asset generator) so PPP and CSF share them; add a CSF dense-window fetcher that reads locus windows by remote tabix from the GWAS Catalog harmonised files; add GRCh38-harmonized trait tasks for Alzheimer's, schizophrenia, Parkinson's and ALS.

**Tech Stack:** Python 3.13, polars, attrs, pysam (TabixFile, remote HTTPS), gwaslab (sumstats + liftover), the repo's Task build system, pytest.

**Spec:** `experiments/claude/design_specs/2026-11-05-csf-protein-coloc-design.md`

## Global Constraints

- All commands run via `pixi r ...`; after each task run `pixi r invoke green 2>&1 | tee <scratch log>` and read the log (testmon may report "no tests ran"; exit 0 alone proves nothing).
- PPP asset ids, metas and on-disk paths unchanged: `{base}_ppp_coloc`, `{base}_ppp_coloc_display`, `{base}_ppp_coloc_plots`, the PPP hit table id/meta (trait `ukbb_ppp`, project `ppp_hit_segments`, sub_dir `analysis`); PPP task `deps` lists unchanged in content **and order**.
- CSF slim asset ids/paths exactly as on `csf-database` (`csf_slim_hapmap_3_<GENE>_<seqid_underscored>`, trait `western_csf_pqtl`, project = gene symbol, `aligned/hapmap_3_index/<asset_id>.parquet`).
- If any already-built asset (CSF slim files, CSF index, PPP coloc outputs, PPP hit table) fails trace verification: **stop and report to the user**; never let it rebuild.
- Protein model for both sources: `QuantitativeTraitModel(sd_y=1.0)`.
- No chrX in CSF.
- Repo conventions: attrs `@frozen`; NewType/Literal ids; no `Any`, no `Callable[...]` (Protocols instead); no monkeypatching — inject; imports at top of file; docstrings without backticks or RST; column-name constants, not literals; read tables via `scan_dataframe_asset`; Task-level tests preferred; no structural tests for asset-wiring modules; never write under `docs/`; never commit logs or PDFs; stage explicit paths.

## Review Focus

1. A CSF window on a chromosome absent from the tabix index (chrX = 23, or any contig not in `TabixFile.contigs`) → empty WINDOW_SCHEMA frame, no exception. Test in Task 6.
2. A remote read failure (unreachable URL, truncated stream) → raises `OSError` so `call_with_retries` retries it. Test in Task 6.
3. `NA` in effect_allele_frequency / beta / standard_error in a harmonised row → null in the frame, not a crash; null-effect rows are later dropped by `_usable_variants`. Test in Task 6.
4. Several aptamers on one gene (APOE has 4) → distinct protein ids end to end, each with its own row in raw results and display. Test in Task 7 (gene coords) and Task 8 (hit table).
5. Multi-gene symbol or symbol missing from MAGMA → no gene interval → cis status "unknown", not a crash, not a wrong "trans". Test in Task 7.

---

### Task 1: Read-only trace verification script

**Files:**
- Create: `experiments/claude/csf_coloc/verify_traces.py`

**Interfaces:**
- Produces: `verify_tasks(tasks: Sequence[Task]) -> VerifyReport` (CLI: `pixi r python experiments/claude/csf_coloc/verify_traces.py <which>` where which ∈ `csf_database`, `ppp_coloc`).

Purpose: tell whether existing assets would verify **without running the runner** (the runner has no dry-run and would start rebuilding). For each task: resolve its path with `DEFAULT_RUNNER.meta_to_path`, check existence, then reproduce `verify_trace` from `mecfs_bio/build_system/rebuilder/verifying_trace_rebuilder/verifying_trace_rebuilder_core.py` read-only: load the info store (`DEFAULT_RUNNER`'s `info_store` path, `VerifyingTraceInfo`), recompute the value trace with `ImoHasher.with_xxhash_128()`, compare with the recorded trace, and compare each recorded dep trace with the dep's current trace (dep resolved by asset id among the transitive deps of the requested tasks). Never write the info store.

- [ ] **Step 1:** Implement with a `VerifyReport` attrs class: counts of `verified`, `missing_file`, `no_trace`, `value_mismatch`, `dep_mismatch`, and up to 20 example asset ids per failure class. Print it.
- [ ] **Step 2:** `ppp_coloc` mode = the four `*_PPP_COLOC` groups' coloc/display/plot tasks + `HAPMAP_3_PPP_HIT_TABLE_5E8`. Run it now, on unmodified main code: `pixi r python experiments/claude/csf_coloc/verify_traces.py ppp_coloc | tee <scratch>/verify_ppp_before.log`. Expected: every task `verified`. Record the result in the commit message; if not all verified, stop and report (the baseline itself is broken).
- [ ] **Step 3:** Commit the script (not the log).

### Task 2: CSF database onto the branch

**Files (from `git show csf-database:<path>`):**
- Create: `mecfs_bio/build_system/task/csf_database/{__init__,build_slim_aptamer_parquet_task,construct_csf_variant_index_task,gwas_catalog_url}.py`
- Create: `mecfs_bio/asset_generator/csf_slim_aptamer_asset_generator.py`
- Create: `mecfs_bio/assets/gwas/csf_pqtl/{__init__,raw/__init__,raw/csf_template_aptamer,csf_database/__init__,csf_database/hapmap3/__init__,csf_database/hapmap3/hapmap3_csf_index,csf_database/hapmap3/hapmap3_csf_database_aptamer_files}.py`
- Create: `mecfs_bio/assets/reference_data/csf_pqtl_sumstats/{__init__.py,csf_aptamer_manifest.csv,regenerate_csf_manifest.py}`
- Create: `mecfs_bio/constants/csf_database_constants.py`, `mecfs_bio/constants/gwas_ssf_constants.py` (if not on main)
- Create: `mecfs_bio/analysis/csf_pqtl_database_build.py`
- Create: `test_mecfs_bio/unit/build_system/task/csf_database/{__init__,test_build_slim_aptamer_parquet_task,test_construct_csf_variant_index_task,test_gwas_catalog_url}.py`
- Modify as needed: `mecfs_bio/build_system/task/ppp_database/allele_key.py` (only what the CSF code imports)
- Do NOT bring: `csf_ldsc/`, `batched_ldsc/`, `hapmap3_csf_heritability.py`, `csf_heritability.py`, CSF LDSC constants/tests/experiments.

**Interfaces:**
- Produces: `HAPMAP_3_CSF_DATABASE: CsfSlimAptamerTaskCollection` (`.aptamer_tasks: tuple[BuildSlimCsfAptamerParquetTask, ...]`, each with `.aptamer: CsfAptamerFile` — fields `analyte, seq_id, accession, uniprot, entrez_gene_symbol`); `HAPMAP_3_CSF_DATABASE_INDEX`; `CSF_APTAMER_MANIFEST_PATH`; `gwas_catalog_bucket(accession) -> str`; `GcstAccession`, `SeqId`, `Analyte`, `UniProtId`; `WESTERN_CSF_TRAIT_NAME = "western_csf_pqtl"`.

- [ ] **Step 1:** Copy the files. Adapt to main: `unordered_allele_key` callers add the `assert_all_snv` guard main now requires; downloads use main's `request_connections` WF/downloader parameter; fix any other import/API drift until `pixi r invoke green` passes. Asset ids, metas, `create()` signatures and path construction unchanged.
- [ ] **Step 2:** Run the brought-over tests: `pixi r pytest test_mecfs_bio/unit/build_system/task/csf_database -p no:testmon -q`. Expected: all pass.
- [ ] **Step 3: USER STEP.** Ask the user to run `pixi r invoke migrate-asset-store --reverse --only western_csf_pqtl` (dry run), then with `--execute`, then delete the `gwas/western_csf_pqtl` line from `default_runner_config.yaml`. Wait for confirmation.
- [ ] **Step 4:** Add `csf_database` mode to the Task 1 script (index task + all 7,008 slim tasks + template aptamer tasks) and run it, tee'd. Expected: 7,008 slim + index `verified`. If not: stop and report (spec rule; a rebuild is ~1.4 TB). In particular check `HAPMAP_3_MEMBERSHIP_LIST`, an upstream dep whose code changed on main.
- [ ] **Step 5:** Commit (code, manifest, tests; no logs).

### Task 3: Neutral protein id and pluggable gene coordinates

**Files:**
- Modify: `mecfs_bio/constants/protein_coloc_constants.py` (add `ProteinId`)
- Create: `mecfs_bio/build_system/task/protein_coloc/gene_coords_reader.py`
- Modify: `mecfs_bio/build_system/task/protein_coloc/candidate_selection.py` (`Oid` → `ProteinId`)
- Modify: `mecfs_bio/build_system/task/protein_coloc/trait_protein_coloc_task.py`
- Modify: `mecfs_bio/asset_generator/ppp_coloc_asset_generator.py` (pass the PPP reader)
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_trait_protein_coloc_task.py`, `test_candidate_selection.py` (update types / construction only)

**Interfaces:**
- Produces:
  - `ProteinId = NewType("ProteinId", str)` in `protein_coloc_constants.py`.
  - In `gene_coords_reader.py`:
    ```python
    class GeneCoordsReader(Protocol):
        def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]: ...

    @frozen(slots=True)
    class St3GeneCoordsReader:  # PPP: wraps ppp_ldsc.gene_coords.read_gene_coords, keys ProteinId(str(oid))
        def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]: ...
    ```
  - `build_candidates(loci, hit_segments, gene_coords: Mapping[ProteinId, GenomicInterval], cis_window_bp)`.
  - `TraitProteinColocTask` gains field `gene_coords_reader: GeneCoordsReader` (after `gene_coords_task`); `create(..., gene_coords_task: Task, gene_coords_reader: GeneCoordsReader, ...)`; execute calls `self.gene_coords_reader(_read_table(self.gene_coords_task, fetch))`. Docstring: gene_coords_task is any table the reader understands. `deps` unchanged.

- [ ] **Step 1:** Update the existing coloc-task tests to construct with `gene_coords_reader=St3GeneCoordsReader()`; run, expect failure (unknown argument).
- [ ] **Step 2:** Implement. `generate_ppp_coloc_tasks` passes `St3GeneCoordsReader()`.
- [ ] **Step 3:** `pixi r pytest test_mecfs_bio/unit/build_system/task/protein_coloc -p no:testmon -q` → all pass, no expected-value changes.
- [ ] **Step 4:** green; commit.

### Task 4: Source-independent hit table

**Files:**
- Modify: `mecfs_bio/build_system/task/protein_coloc/protein_hit_table_task.py`
- Modify: `mecfs_bio/assets/gwas/ukbb_ppp/ppp_database/hapmap3/eur_discovery_hapmap3_ppp_hit_table.py`
- Modify: `mecfs_bio/build_system/task/protein_coloc/trait_protein_coloc_task.py` (field type)
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_protein_hit_table_task.py`

**Interfaces:**
- Produces:
  ```python
  @frozen(slots=True)
  class HitTableProtein:
      task: Task          # slim per-protein BETA/SE file aligned to the index
      protein_id: ProteinId
      gene: str

  @frozen(slots=True)
  class ProteinHitTableTask(GeneratingTask):
      meta: Meta
      proteins: tuple[HitTableProtein, ...]
      index_task: Task
      p_threshold: float
      merge_gap_bp: int
      # deps == [index_task, *(p.task for p in proteins)]  (same order as before)

      @classmethod
      def create(cls, asset_id: str, proteins: tuple[HitTableProtein, ...], index_task: Task,
                 trait: str, project: str, p_threshold: float = 5e-8,
                 merge_gap_bp: int = 250_000) -> "ProteinHitTableTask": ...

  def ppp_hit_table_proteins(tasks: Sequence[BuildSlimProteinParquetTask]) -> tuple[HitTableProtein, ...]
      # ProteinId(str(t.protein.oid)), gene=t.protein.gene; lives in the PPP hit-table asset module
  ```
  `PppProteinHitTableTask` is removed (renamed); `TraitProteinColocTask.hit_table_task: ProteinHitTableTask`.
- PPP asset passes `trait="ukbb_ppp", project="ppp_hit_segments"`; meta otherwise identical to today (`ResultTableMeta`, sub_dir `analysis`, `.parquet`, parquet read spec).

- [ ] **Step 1:** Update hit-table tests to the new constructor; add `test_hit_table_keeps_aptamers_of_one_gene_distinct`: two `HitTableProtein`s with the same `gene="APOE"` and protein ids `"2418-55"` / `"2937-10"`, each with a significant variant → table has two rows with those two `HIT_PROTEIN_ID_COL` values. Run, expect failure.
- [ ] **Step 2:** Implement; update the PPP asset module.
- [ ] **Step 3:** Tests pass (`... protein_coloc -p no:testmon -q`).
- [ ] **Step 4:** green; commit.

### Task 5: Shared coloc asset generator; prove PPP unchanged

**Files:**
- Create: `mecfs_bio/asset_generator/protein_coloc_asset_generator.py`
- Modify: `mecfs_bio/asset_generator/ppp_coloc_asset_generator.py` (thin wrapper)
- Modify: `mecfs_bio/analysis/ppp_coloc_analysis.py` (type name only)

**Interfaces:**
- Produces:
  ```python
  @frozen(slots=True)
  class ProteomicsColocSource:
      asset_infix: Literal["ppp", "csf"]
      hit_table_task: ProteinHitTableTask
      gene_coords_task: Task
      gene_coords_reader: GeneCoordsReader
      fetcher: DenseWindowFetcher

  @frozen(slots=True)
  class ColocTaskGroup:  # renamed from PppColocTaskGroup
      coloc_task: TraitProteinColocTask
      display_task: ColocDisplayTableTask
      plot_task: ColocLocusPlotTask
      def terminal_tasks(self) -> list[Task]: ...

  def generate_protein_coloc_tasks(base_name: str, source: ProteomicsColocSource,
      trait_task: Task, trait_model: TraitModel, config: ColocRunConfig,
      h4_threshold: float = 0.8) -> ColocTaskGroup
      # asset ids f"{base_name}_{infix}_coloc", f"..._coloc_display", f"..._coloc_plots";
      # FASTA UCSC_HG38_INDEXED_FASTA, panel GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES "nfe",
      # gene track MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED

  def generate_ppp_coloc_tasks(base_name, trait_task, trait_model, config=ColocRunConfig(),
      h4_threshold=0.8, hit_table_task=HAPMAP_3_PPP_HIT_TABLE_5E8,
      manifest_path=EUR_DISCOVERY_PPP_MANIFEST_PATH) -> ColocTaskGroup  # signature unchanged
  ```
- [ ] **Step 1:** Implement.
- [ ] **Step 2:** Run `verify_traces.py ppp_coloc` again, tee'd. Expected: identical to Task 1 Step 2 (all verified). Also compare asset ids: print the 12 PPP coloc asset ids + hit table id before (from Task 1 log) and after — identical. If anything differs: stop and report.
- [ ] **Step 3:** green; commit.

### Task 6: CSF dense-window fetcher (remote tabix)

**Files:**
- Modify: `mecfs_bio/build_system/task/csf_database/gwas_catalog_url.py` (add harmonised URL)
- Create: `mecfs_bio/build_system/task/protein_coloc/csf_window_fetcher.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_csf_window_fetcher.py`, extend `test_gwas_catalog_url.py`

**Interfaces:**
- Produces:
  ```python
  def gwas_catalog_harmonised_url(accession: GcstAccession) -> str
      # f"{_FTP_ROOT}/{bucket}/{accession}/harmonised/{accession}.h.tsv.gz"

  class HarmonisedLocationResolver(Protocol):
      def __call__(self, accession: GcstAccession) -> str: ...   # URL or local path

  @frozen(slots=True)
  class GwasCatalogHarmonisedResolver:  # production default -> gwas_catalog_harmonised_url
      def __call__(self, accession: GcstAccession) -> str: ...

  @frozen(slots=True)
  class CsfWindowFetcher:  # implements DenseWindowFetcher
      aptamers: tuple[CsfAptamerFile, ...]     # looked up by ProteinId(seq_id)
      resolver: HarmonisedLocationResolver = GwasCatalogHarmonisedResolver()
      def fetch(self, protein_id: str, windows: Sequence[LocusWindow],
                scratch_dir: Path, wf: WF) -> dict[str, pl.DataFrame]: ...
  ```
- Column contract (verified 2026-11-05, constant tuple in the module, assert each row has 15 fields): `chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value variant_id rsid neg_log_10_p_value rs_id n hm_coordinate_conversion hm_code`. Map to WINDOW_SCHEMA: CHR←chromosome, POS←base_pair_location, EA←effect_allele, NEA←other_allele, BETA←beta, SE←standard_error, EAF←effect_allele_frequency, N←n; `NA` → null. Contig name is `str(window.chrom)`; contigs not in `TabixFile.contigs` → empty frame. pysam fetch bounds are 0-based half-open: `fetch(contig, start - 1, end)` for the closed 1-based window. Wrap pysam open/fetch errors so failures surface as `OSError` (pysam raises `OSError`/`ValueError`; re-raise non-OSError I/O failures as `OSError`). Downloads nothing to scratch.
- Tests build a small bgzipped, tabix-indexed fixture in `tmp_path` with `pysam.tabix_compress` + `pysam.tabix_index(..., seq_col=0, start_col=1, end_col=1, line_skip=1)` and a header line, injected via a test resolver returning the local path.

- [ ] **Step 1: Failing tests:**
  - `test_fetch_returns_window_rows_in_window_schema`: fixture rows on chr 19 at 100, 200, 300; window [150, 300] → POS [200, 300]; columns/dtypes == WINDOW_SCHEMA; N values from the `n` column.
  - `test_fetch_window_on_absent_contig_is_empty`: window on chrom 23 → empty frame with WINDOW_SCHEMA (Review Focus 1).
  - `test_fetch_na_fields_become_null`: a row with `NA` EAF and `NA` beta → nulls (Review Focus 3).
  - `test_fetch_unreachable_source_raises_oserror`: resolver returns a nonexistent path → `pytest.raises(OSError)` (Review Focus 2).
  - `test_fetch_unknown_protein_id_asserts`: protein id not in `aptamers` → `AssertionError`.
  - `test_gwas_catalog_harmonised_url` for GCST90421540 → `.../GCST90421001-GCST90422000/GCST90421540/harmonised/GCST90421540.h.tsv.gz`.
- [ ] **Step 2:** Run, expect failures.
- [ ] **Step 3:** Implement.
- [ ] **Step 4:** Tests pass.
- [ ] **Step 5:** Smoke (not a unit test, experiments script `experiments/claude/csf_coloc/fetcher_smoke.py`): real fetch for APOE aptamer seq id 2418-55 window chr19:44,400,000-45,400,000 → a few thousand rows. Tee the log.
- [ ] **Step 6:** green; commit.

### Task 7: CSF gene coordinates from the manifest and MAGMA

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/csf_gene_coords_reader.py`
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_csf_gene_coords_reader.py`

**Interfaces:**
- Produces:
  ```python
  @frozen(slots=True)
  class CsfGeneCoordsReader:  # implements GeneCoordsReader
      aptamers: tuple[CsfAptamerFile, ...]
      def __call__(self, table: pl.DataFrame) -> dict[ProteinId, GenomicInterval]: ...
  ```
  `table` is `MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED` read through its meta (GENE_INFO_NAME_COL, GENE_INFO_CHROM_COL, GENE_INFO_START_COL, GENE_INFO_END_COL; chromosome "X" → 23 via `parse_chrom`). Keys `ProteinId(seq_id)`. A symbol containing a multi-gene separator (whatever the manifest uses for its 78 multi-gene aptamers — inspect them), or absent from MAGMA, or present on more than one MAGMA row, → no entry.
- [ ] **Step 1: Failing tests:**
  - `test_symbol_found_gives_interval`
  - `test_two_aptamers_same_gene_both_mapped` (Review Focus 4): seq ids 2418-55 and 2937-10, both APOE → two keys, same interval.
  - `test_multi_gene_symbol_has_no_interval` (Review Focus 5)
  - `test_symbol_missing_from_magma_has_no_interval` (Review Focus 5)
  - `test_duplicate_magma_symbol_has_no_interval`
  - A coloc-task-level test in `test_trait_protein_coloc_task.py`, `test_csf_unmapped_aptamer_is_cis_unknown`: run the task with a CSF-shaped source whose aptamer symbol is missing from the gene table → raw results `CIS_STATUS_COL == CIS_STATUS_UNKNOWN` for that pair (reuse the existing test's fake fetcher pattern).
- [ ] **Step 2:** Run, expect failures. **Step 3:** Implement. **Step 4:** Pass. **Step 5:** green; commit.

### Task 8: CSF hit table asset and CSF coloc generator

**Files:**
- Create: `mecfs_bio/assets/gwas/csf_pqtl/csf_database/hapmap3/hapmap3_csf_hit_table.py`
- Create: `mecfs_bio/asset_generator/csf_coloc_asset_generator.py`
- Test: extend `test_protein_hit_table_task.py` only if Task 4's test does not already cover CSF-shaped input.

**Interfaces:**
- Produces:
  - `HAPMAP_3_CSF_HIT_TABLE_5E8 = ProteinHitTableTask.create(asset_id="hapmap3_csf_hit_table_5e8", proteins=csf_hit_table_proteins(HAPMAP_3_CSF_DATABASE.aptamer_tasks), index_task=HAPMAP_3_CSF_DATABASE_INDEX, trait="western_csf_pqtl", project="csf_hit_segments")`; `csf_hit_table_proteins` maps each task to `HitTableProtein(task, ProteinId(seq_id), gene=entrez_gene_symbol)`.
  - `CSF_SELECTION = SelectionConfig(max_pairs=200, max_per_locus=20)`; `CSF_COLOC_CONFIG = ColocRunConfig(selection=CSF_SELECTION)`. Docstring states why: one CSF pair costs one ~4 s remote tabix read (200 pairs ≈ 15 min) versus a Synapse tar per PPP protein; per-locus cap raised because several aptamers can share one gene.
  - `generate_csf_coloc_tasks(base_name: str, trait_task: Task, trait_model: TraitModel, config: ColocRunConfig = CSF_COLOC_CONFIG, h4_threshold: float = 0.8) -> ColocTaskGroup` using `ProteomicsColocSource(asset_infix="csf", hit_table_task=HAPMAP_3_CSF_HIT_TABLE_5E8, gene_coords_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED, gene_coords_reader=CsfGeneCoordsReader(aptamers), fetcher=CsfWindowFetcher(aptamers))`, aptamers from the manifest.
- Display / plot: the GENE column comes from the hit table's gene; PROTEIN_ID is the seq id, so aptamers of one gene stay distinguishable. No code change expected; confirm with the Task 7 task-level test output.
- [ ] **Step 1:** Implement both modules.
- [ ] **Step 2:** Materialize the hit table: `DEFAULT_RUNNER.run([HAPMAP_3_CSF_HIT_TABLE_5E8], ...)` from a short script in `experiments/claude/csf_coloc/`, tee'd. Expected: slim files verify (no downloads logged), table written; log rows, distinct protein ids, and the APOE rows.
- [ ] **Step 3:** green; commit.

### Task 9: Alzheimer's and Parkinson's GRCh38 harmonized traits (no liftover)

**Files:**
- Create: `mecfs_bio/assets/gwas/alzheimers/bellenguez_et_al/processed/bellenguez_build_38_harmonized.py`
- Create: `mecfs_bio/assets/gwas/parkinsons/nalls_et_al_2019/{__init__,raw/__init__,raw/nalls_2019_raw,processed/__init__,processed/nalls_2019_build_38_harmonized}.py`

Pattern: `willer_ldl_build_38_harmonized.py` (GWASLabCreateSumstatsTask genome_build="38" → GwasLabSumstatsToTableTask → GenomeReferenceHarmonizationTask with UCSC_HG38_INDEXED_FASTA, GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES, "nfe").

- **AD** from `BELLENGUEZ_ET_AL_ALZHEIMERS_RAW` (GCST90027158 harmonised): chrom `hm_chrom`, pos `hm_pos`, ea `hm_effect_allele`, nea `hm_other_allele`, beta `hm_beta`, se `standard_error`, p `p_value`, eaf `hm_effect_allele_frequency`, n = `n_cas + n_con` (ExprPipe into GWASLAB_SAMPLE_SIZE_COLUMN); DropNullsPipe on hm_chrom/hm_pos. Asset id `bellenguez_alzheimers_build_38_genome_reference_harmonized`.
- **PD** download `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST009001-GCST010000/GCST009325/harmonised/GCST009325.h.tsv.gz`, md5 `e7f2654a34931bb9996fb0927cece282`. This new-format harmonised file is GRCh38 (`hm_coordinate_conversion`); columns `chromosome base_pair_location effect_allele other_allele beta standard_error effect_allele_frequency p_value rsid N_cases N_controls hm_coordinate_conversion hm_code variant_id`; n = `N_cases + N_controls`. Asset id `nalls_2019_parkinsons_build_38_genome_reference_harmonized`. Docstring: Nalls et al. 2019 without 23andMe, 15,056 cases + 18,618 UK Biobank proxy cases + 449,056 controls, EUR.
- [ ] **Step 1:** Implement both modules (no unit tests — asset wiring).
- [ ] **Step 2:** Materialize both harmonized tasks (script in `experiments/claude/csf_coloc/`, tee'd). Check: PD position of rs356182 equals its GRCh38 position per Ensembl (GRCh38 chr4:89,704,960; GRCh37 would be 90,626,111) — confirms the file is GRCh38; AD has genome-wide-significant variants near APOE (chr19 ~44.9 Mb). If the PD file is not GRCh38, stop and switch PD to the raw-file liftover route of Task 10.
- [ ] **Step 3:** green; commit.

### Task 10: Schizophrenia and ALS GRCh38 harmonized traits (gwaslab liftover 37 → 38)

**Files:**
- Create: `mecfs_bio/assets/gwas/schizophrenia/pgc2022/processed/pgc_2022_sch_build_38_harmonized.py`
- Create: `mecfs_bio/assets/gwas/als/van_rheenen_et_al_2021/{__init__,raw/__init__,raw/van_rheenen_2021_raw,processed/__init__,processed/van_rheenen_2021_build_38_harmonized}.py`

Pattern: `GWASLabCreateSumstatsTask(genome_build="19", liftover_to="38", ...)` (liftover usage exists in `liu_et_al_2023_eur_liftover_to_37_sumstats.py`, other direction) → table → `GenomeReferenceHarmonizationTask` as in Task 9.

- **SCZ** from `PGC_2022_SCH_RAW` (GRCh37; the user believes no native GRCh38 release exists): chrom `CHROM`, pos `POS`, ea `A1`, nea `A2`, beta `BETA` (log OR), se `SE`, p `PVAL`, eaf `FCON`, n = `NCAS + NCON`. Asset id `pgc_2022_sch_build_38_genome_reference_harmonized`.
- **ALS**: the harmonised file lacks `N_effective`, so use the raw file: `https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90027001-GCST90028000/GCST90027164/GCST90027164_buildGRCh37.tsv.gz` (md5 from that directory's `md5sum.txt`; fetch it while implementing and pin it). Columns `rsid chromosome base_pair_location effect_allele other_allele effect_allele_frequency beta standard_error p_value N_effective`; alleles lowercase → uppercase (ExprPipe as Willer); n = `N_effective`. Asset id `van_rheenen_2021_als_eur_build_38_genome_reference_harmonized`.
- Memory: 16 GB machine; if gwaslab OOMs, add a MAF ≥ 0.01 filter first (RA precedent) and say so in the docstring.
- [ ] **Step 1:** Implement.
- [ ] **Step 2:** Materialize both (tee'd). Check lifted positions for one known SNP each against Ensembl GRCh38, and that GWS loci exist (SCZ: MHC region chr6 and many others; ALS: C9orf72 region chr9 ~27.5 Mb, UNC13A chr19 ~17.6 Mb).
- [ ] **Step 3:** green; commit.

### Task 11: Per-trait CSF coloc modules and driver

**Files:**
- Create: `.../analysis/<trait>_csf_coloc.py` next to each existing `<trait>_ppp_coloc.py` (DecodeME, RA, LDL, IBD); new `analysis/` modules for AD (`bellenguez_alzheimers_csf_coloc.py`), SCZ (`pgc_2022_sch_csf_coloc.py`), PD (`nalls_2019_parkinsons_csf_coloc.py`), ALS (`van_rheenen_2021_als_csf_coloc.py`).
- Create: `mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/decode_me_csf_coloc_exploratory.py`
- Create: `mecfs_bio/analysis/csf_coloc_analysis.py`

**Trait models (each via `generate_csf_coloc_tasks`):**
- DecodeME, RA, LDL, IBD: exactly the trait_model of their `*_ppp_coloc.py`; base_name the same as their PPP module.
- AD: `CaseControlTraitModel(case_fraction=0.1763)` — (39,106 clinical + 46,828 proxy cases) / 487,511 (GWAS Catalog GCST90027158); proxies counted as cases.
- SCZ: existing `prevalance_info.py` `sample_prevalence` (0.408).
- PD: `CaseControlTraitModel(case_fraction=0.0698)` — (15,056 + 18,618 proxy) / 482,730; proxies counted as cases (spec).
- ALS: `CaseControlTraitModel(case_fraction=0.5)` with N = N_effective (spec).
- Exploratory DecodeME: two groups, base names `decode_me_gwas_1_lead_p_1e_6` and `decode_me_gwas_1_lead_p_1e_5`, `ColocRunConfig(lead_p_threshold=1e-6 | 1e-5, selection=SelectionConfig(max_pairs=1000, max_per_locus=20))`, docstring "exploratory: sub-genome-wide-significant DecodeME loci". Before fixing max_pairs, count DecodeME loci at each threshold (`find_trait_loci` on the harmonized trait, small experiments script, tee'd) and raise max_pairs if 1000 would leave candidates unselected at under ~3 per locus; record the counts in the docstring.
- Driver: `run_<trait>_csf_coloc()` per group plus `run_decode_me_csf_coloc_exploratory()`; `main()` runs AD first only (the gate), others are called explicitly after the gate.
- [ ] **Step 1:** Implement. **Step 2:** green (imports construct every task). **Step 3:** commit.

### Task 12: Alzheimer's positive-control gate

**Files:**
- Create: `experiments/claude/csf_coloc/western_ad_proteins.csv` (+ a provenance header comment or sibling `.md`)
- Create: `experiments/claude/csf_coloc/ad_positive_control.py`

- [ ] **Step 1:** Obtain the Western et al. 2024 paper and supplementary tables (Nature Genetics 56:2672, PMID 39528825). Read the PDF/tables (`pixi r -e pdf-env pdftotext`, openpyxl for .xlsx) — do not work from the abstract or memory. Extract the proteins/aptamers they report as colocalizing with / causal for AD (seq id, gene, their PP.H4 or MR evidence, table and page). Never commit the PDF.
- [ ] **Step 2:** Run the AD group (`run_bellenguez_alzheimers_csf_coloc`), tee'd.
- [ ] **Step 3:** `ad_positive_control.py` joins the reference list with our raw results on seq id and reports per protein: candidate (hit segment overlapped an AD locus)? selected within budget? tested? PP.H4; and lists our PP.H4 ≥ 0.8 pairs not in their list. Tee to log.
- [ ] **Step 4:** **USER GATE.** Report the comparison to the user. Misses not explained by method differences or budget are investigated (fetcher, protein ids, orientation) before Task 13. The user judges pass/fail.
- [ ] **Step 5:** Commit the reference CSV and script (not logs).

### Task 13: Remaining runs

- [ ] **Step 1:** After the user passes the gate: run the seven other trait groups and the two exploratory DecodeME groups (tee'd logs).
- [ ] **Step 2:** Report per trait: loci, candidates, tested, failed, PP.H4 ≥ 0.8 pairs (gene, seq id, cis/trans).
