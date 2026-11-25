# Protein Coloc Follow-ups Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:**
- Cross-check our coloc.abf against R coloc.
- Add context columns to the coloc display table.
- Add per-locus stacked Manhattan plots and a protein-annotation table.
- Run all of it on the exploratory DecodeME x PPP coloc run.

**Architecture:**
- Everything derives from an existing TraitProteinColocTask result directory: raw_results.parquet plus windows/*.parquet holding the aligned betas. No coloc reruns.
- Two new Tasks (plots and info) select loci from the display table with one shared ranking helper.
- Plot helpers are hoisted into a shared module. The R check is an experiments script running R coloc in a throwaway `pixi exec` environment.

**Tech Stack:**
- polars and matplotlib (Figure API).
- gget, and the Ensembl REST API via requests.
- R coloc 5.2.3, via `pixi exec`.
- attrs @frozen Tasks.

**Spec:** experiments/claude/design_specs/2026-11-13-coloc-followups-design.md

## Global Constraints

- Run everything via `pixi r ...`.
- After significant changes run `pixi r invoke green 2>&1 | tee <scratchpad>/green_<task>.log`, and check the log, not just the exit code.
- Never commit logs or PDFs. Stage explicit paths only. Never edit docs/.
- No PPP or CSF rebuilds, and TraitProteinColocTask is not re-run.
- Display column names stay short: `trait_p`, `prot_p`, `same_sign`, `n_snp`.
- Docstrings: no backticks or RST. Define acronyms.
- Imports go at the top of the file.
- Use column constants, not repeated literals.
- Inject dependencies instead of monkeypatching.
- Prefer Task-level tests.
- Commits are local, on branch decode_me_coloc. Do not push.

## Review Focus

- **Empty display table (no PP.H4 > threshold).** Plot and info tasks produce an empty index or table, not a crash. This is tested in Tasks 4 and 5.
- **Multi-gene protein symbol (e.g. "IL12A_IL12B").** The info task splits it into one row per gene. This is tested in Task 5.
- **Symbol that Ensembl cannot resolve.** The row keeps null annotation columns and a warning is logged. This is tested in Task 5.
- **Protein window whose lead beta sign differs from the trait's.** `same_sign` is False. This is tested in Task 2.
- **Locus window extending past the genetic-map coverage (e.g. a chr start at 0).** The recombination track draws whatever rows exist and does not fail on an empty slice. The existing `_plot_recomb` already guards `len(bp) >= 1`; Task 3 keeps that guard.

---

### Task 1: R coloc cross-check script

**Files:**
- Create: `experiments/claude/protein_coloc/r_coloc_crosscheck.py`
- Create: `experiments/claude/protein_coloc/r_coloc_crosscheck.R`

**Interfaces:**
- Consumes:
  - `DECODE_ME_PPP_COLOC_EXPLORATORY_2.coloc_task` and its `.trait_model.case_fraction`.
  - The coloc result dir, located via the task meta (`simple_meta_to_relative_path` plus the runner config). Do not use `find`.
  - Constants: `RAW_RESULTS_FILENAME`, `WINDOWS_DIRNAME`, `pair_file_stem`, `PP_H0_COL`..`PP_H4_COL`, `P1_COL`, `P2_COL`, `P12_COL`, `STATUS_TESTED`.
- Produces: a log at `experiments/claude/logs/r_coloc_crosscheck.log` (uncommitted) and a summary section appended to `experiments/claude/protein_coloc/v1_run_summary.md`.

- [ ] **Step 1: Write the R script.**
  - Arguments: an input dir of per-pair CSVs plus a priors CSV, and an output CSV path.
  - For each pair it calls:

    ```r
    coloc::coloc.abf(
      dataset1 = list(beta=d$BETA_TRAIT, varbeta=d$SE_TRAIT^2, type="cc", s=case_fraction, snp=d$snp),
      dataset2 = list(beta=d$BETA_PROTEIN, varbeta=d$SE_PROTEIN^2, type="quant", sdY=1, snp=d$snp),
      p1=p1, p2=p2, p12=p12)$summary
    ```

  - It writes PP.H0..PP.H4 per pair.
  - `snp` is `paste(CHR, POS, EA, NEA)`.
  - Before writing the comparison, confirm by reading coloc 5.2.3's `process.dataset`/`approx.bf.estimates` source that:
    - for type cc with beta and varbeta given, sd.prior is 0.2 and MAF is not required;
    - for quant with sdY = 1, sd.prior is 0.15.
  - Note the answer in a script comment.

- [ ] **Step 2: Write the Python driver.**
  - It exports every tested pair's window to CSV in a temp dir under the scratchpad.
  - It runs:

    ```
    pixi exec --spec r-base=4.4 --spec r-coloc=5.2.3 Rscript r_coloc_crosscheck.R <in> <out>
    ```

    through `execute_command`.
  - It joins the result to raw_results on (LOCUS_ID, PROTEIN_ID).
  - It prints:
    - max |dPP| per hypothesis over all tested pairs;
    - the same for the display-table hits (PP_H4 > 0.8);
    - the 5 worst pairs.
  - Exit non-zero if max |dPP| >= 1e-6.

- [ ] **Step 3: Run.**
  - Command: `pixi r python experiments/claude/protein_coloc/r_coloc_crosscheck.py 2>&1 | tee experiments/claude/logs/r_coloc_crosscheck.log`
  - Expected: 517 pairs compared, max |dPP| < 1e-6.
  - If this fails, stop and report to the user with the worst pairs. Do not proceed to "fix" either side without discussion.

- [ ] **Step 4: Append the result summary to `v1_run_summary.md` and commit.**

  ```bash
  git add experiments/claude/protein_coloc/r_coloc_crosscheck.py experiments/claude/protein_coloc/r_coloc_crosscheck.R experiments/claude/protein_coloc/v1_run_summary.md
  git commit -m "Cross-check Python coloc.abf against R coloc 5.2.3 on DecodeME exploratory run"
  ```

### Task 2: Display table context columns

**Files:**
- Modify: `mecfs_bio/constants/protein_coloc_constants.py`. Add `DISPLAY_TRAIT_P_COL = "trait_p"`, `DISPLAY_PROTEIN_P_COL = "prot_p"`, `DISPLAY_SAME_SIGN_COL = "same_sign"` and `DISPLAY_N_SNP_COL = "n_snp"`.
- Modify: `mecfs_bio/build_system/task/protein_coloc/coloc_display_table_task.py`.
- Modify: `test_mecfs_bio/unit/build_system/task/protein_coloc/coloc_result_fixtures.py`, if its windows need controllable betas.
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_display_table_task.py`.

**Interfaces:**
- Produces:
  - `DISPLAY_SCHEMA` gains, after `fdr_sig`: `trait_p: Float64`, `prot_p: Float64`, `same_sign: Boolean`, `n_snp: Int32`.
  - `display_rows(raw: pl.DataFrame, windows_dir: Path, h4_threshold: float, fdr_alpha: float) -> pl.DataFrame`.
  - `ranked_display_loci(display: pl.DataFrame, max_loci: int) -> list[str]`: distinct `locus` values ordered by their max `pp_h4` descending (ties broken by locus string), truncated to max_loci. Tasks 4 and 5 use it.

- [ ] **Step 1: Write the failing tests.**
  - `test_display_adds_lead_p_values_and_n_snp`:
    - set the raw row with `LEAD_MLOG10P=7.0`, `PROTEIN_LEAD_MLOG10P=10.0`, `N_ALIGNED=1234`;
    - assert `trait_p ≈ 1e-7`, `prot_p ≈ 1e-10`, `n_snp == 1234`.
  - `test_same_sign_true_and_false`: two hit pairs whose windows have protein-lead rows (the max |BETA_PROTEIN/SE_PROTEIN|):
    - in one pair BETA_TRAIT and BETA_PROTEIN are both positive;
    - in the other BETA_TRAIT is negative and BETA_PROTEIN positive;
    - assert `same_sign == [True, False]` in PP.H4 order.
  - `test_ranked_display_loci_orders_by_best_pp_h4_and_truncates`:
    - loci "A" (best 0.9), "B" (best 0.95) and "C" (0.85);
    - `ranked_display_loci(df, 2) == ["B", "A"]`.

- [ ] **Step 2: Run them and verify they fail.**
  - Command: `pixi r pytest test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_display_table_task.py -v -p no:testmon`

- [ ] **Step 3: Implement.**
  - p-values: `10 ** -mlog10p` expressions on `LOCUS_LEAD_MLOG10P_COL` and `PROTEIN_LEAD_MLOG10P_COL`.
  - `same_sign`: a helper `_lead_signs_agree(windows_dir: Path, locus_id: str, protein_id: str) -> bool`. It reads `windows_dir / f"{pair_file_stem(...)}.parquet"`, takes the row with max |BETA_PROTEIN/SE_PROTEIN|, and compares `sign(BETA_TRAIT) == sign(BETA_PROTEIN)`.
    - Apply it only to rows surviving the h4 filter.
    - Build the column with a list comprehension in a helper, not in execute.
  - `execute`:
    - fetch the coloc asset once (it is a DirectoryAsset);
    - pass `asset.path / WINDOWS_DIRNAME`;
    - keep `read_raw_results` for the plot task's use.
  - Add the four columns to the module docstring.

- [ ] **Step 4: Run the tests and verify they pass.**
  - Run the test file again.
  - Then run `pixi r invoke green 2>&1 | tee …/green_display.log` and expect 0 failures in the log.

- [ ] **Step 5: Commit** the constants, the task, the tests and the fixtures: `"Add lead p-values, effect-sign agreement and variant count to coloc display table"`.

### Task 3: hg38 genetic map and shared locus-plot helpers

**Files:**
- Create: `mecfs_bio/assets/reference_data/genetic_map/genetic_map_hg38.py`. Mirror `genetic_map_hg19.py`:
  - the URL is `https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg38_withX.txt.gz`;
  - the ids are `genetic_map_hg38_eagle` and `genetic_map_hg38`;
  - `md5_hash` is filled after the first download (Step 3).
- Create: `mecfs_bio/build_system/task/locus_plot_helpers.py`.
- Modify: `mecfs_bio/build_system/task/polyfun_explain/polyfun_explain_plot_task.py` to import `load_recomb` and `plot_recomb` from the new module.
- Modify: `mecfs_bio/build_system/task/protein_coloc/coloc_locus_plot_task.py` to import `plot_association` and `plot_gene_track_panel` from the new module.

**Interfaces:**
- Produces, in `locus_plot_helpers.py`:
  - `load_recomb(fetch: Fetch, task: Task, chrom: int, bp_min: int, bp_max: int) -> pl.DataFrame`: unchanged behaviour. The docstring says "from a genetic map" with no build named.
  - `plot_recomb(ax: Axes, recomb: pl.DataFrame) -> None`: unchanged.
  - `plot_association(ax: Axes, pos: np.ndarray, mlog10p: np.ndarray, label: str, significance_line: bool) -> int`. It scatters, marks and returns the lead index, and draws the dashed `GWAS_SIGNIFICANCE_MLOG10P` line only when `significance_line` is set.
  - `plot_gene_track_panel(ax: Axes, genes: pl.DataFrame, chrom: int, start: int, end: int) -> None`: the body of the current `_plot_genes`.
  - `GENETIC_MAP_HG38: PipeDataFrameTask`, in the new asset file.

- [ ] **Step 1: Move the helpers and update both call sites.**
  - `ColocLocusPlotTask` keeps `significance_line=True` on both of its panels, so its output is unchanged.

- [ ] **Step 2: Run the existing tests for both tasks and verify they pass.**
  - Command: `pixi r pytest test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_locus_plot_task.py test_mecfs_bio/unit/build_system/task/polyfun_explain -v -p no:testmon`
  - Expected: PASS. This is a pure refactor.

- [ ] **Step 3: Materialize GENETIC_MAP_HG38 once and record its md5 in the asset.**
  - Get the md5 with `md5sum` on the downloaded file, located via its meta.
  - Check that the parquet has CHR in 1..23 and POS sorted.

- [ ] **Step 4: Run green (tee to a log) and commit:** `"Add hg38 Eagle genetic map; share locus plot helpers"`.

### Task 4: ColocLocusStackPlotTask

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/coloc_locus_stack_plot_task.py`.
- Modify: `mecfs_bio/build_system/task/protein_coloc/trait_protein_coloc_task.py`. Rename `_scan_trait` to the public `scan_trait`, and update its call site.
- Modify: `mecfs_bio/constants/protein_coloc_constants.py`. Add `N_PROTEINS_COL = "N_PROTEINS"`.
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_locus_stack_plot_task.py`.

**Interfaces:**
- Consumes:
  - `ranked_display_loci` and the display constants (Task 2);
  - `load_recomb`, `plot_recomb`, `plot_association` and `plot_gene_track_panel` (Task 3);
  - `scan_trait(trait_task: Task, fetch: Fetch) -> pl.LazyFrame` (WINDOW_SCHEMA);
  - `read_raw_results`.
- Produces:
  - `ColocLocusStackPlotTask` with these fields:
    - `meta: ResultDirectoryMeta`;
    - `coloc_task: TraitProteinColocTask`;
    - `display_task: ColocDisplayTableTask`;
    - `gene_info_task: Task`;
    - `genetic_map_task: Task`;
    - `max_loci: int`.
  - `create(asset_id, coloc_task, display_task, gene_info_task, genetic_map_task, max_loci=10)`. Meta trait and project come from `coloc_task.meta`. Assert `max_loci >= 1`.
  - Output dir: `<stem>.png` and `<stem>.svg` per locus, where the stem is LOCUS_ID with ":" replaced by "_". Plus `plot_index.parquet` with columns `LOCUS_ID: String`, `N_PROTEINS: Int32`, `PNG_FILENAME: String` and `SVG_FILENAME: String`.

- [ ] **Step 1: Write the failing tests.** The fixtures are:
  - a synthetic coloc dir (reuse `write_coloc_dir`) with 3 loci of hits;
  - a display parquet written via `ColocDisplayTableTask`'s `display_rows`, or by hand with `DISPLAY_SCHEMA`;
  - a tiny trait parquet whose FakeTask meta is readable by `scan_dataframe_asset`;
  - a gene table and a genetic-map parquet.

  Tests:
  - `test_one_figure_per_display_locus`: with `max_loci=10`, the index has 3 rows, each png/svg file exists, and `N_PROTEINS` matches the display rows per locus.
  - `test_max_loci_truncates_by_best_pp_h4`: with `max_loci=1`, the index's single LOCUS_ID is the locus with the highest pp_h4.
  - `test_empty_display_table_gives_empty_index`: the index has 0 rows and no png files.

- [ ] **Step 2: Run them and verify they fail** (import error).

- [ ] **Step 3: Implement.**
  - `execute` reads: display, then `ranked_display_loci`, then per-locus specs from raw_results, then `_render_locus` for each locus, then the index.
  - Loops and file I/O live in helpers.
  - Figure layout:
    - height is `_PANEL_HEIGHT_IN * (2 + n_proteins)`, with `sharex`;
    - the x range is the locus [START, END].
  - Trait row:
    - the trait window comes from `scan_trait(...).filter(CHR == chrom, POS in [start, end])`;
    - -log10 p is from `two_sided_mlog10p(beta/se)`;
    - `significance_line=True`;
    - the recombination rate goes on `ax.twinx()`.
  - Protein rows:
    - data is the saved window for (locus, protein) in display (pp_h4) order;
    - `significance_line=False`;
    - add `ax.axvline(trait_lead_pos, linestyle=":", color="grey")`;
    - the y-label is the gene; the panel title is `f"{gene} ({cis}) PP.H4={pp_h4:.2f}"`.
  - The last row is the gene track panel.
  - The suptitle is the locus id.

- [ ] **Step 4: Run the tests and verify they pass.** Then run green (tee to a log).

- [ ] **Step 5: Commit:** `"Add per-locus stacked trait/protein Manhattan plot task for coloc results"`.

### Task 5: ColocProteinInfoTask

**Files:**
- Create: `mecfs_bio/build_system/task/protein_coloc/coloc_protein_info_task.py`.
- Modify: `mecfs_bio/build_system/task/fetch_gget_info_task.py`. Extract `fetch_gget_info(ensembl_ids: list[str]) -> pd.DataFrame`, which is the retry plus result normalisation currently inline in `execute`, and make `FetchGGetInfoTask` call it.
- Test: `test_mecfs_bio/unit/build_system/task/protein_coloc/test_coloc_protein_info_task.py`.

**Interfaces:**
- Consumes: `ranked_display_loci`, the display constants, and `GGET_TRANSIENT_ERRORS`.
- Produces:
  - `SymbolLookup` (Protocol): `__call__(self, symbol: str) -> str | None`. It returns an Ensembl gene id, or None when unresolved.
  - `InfoLookup` (Protocol): `__call__(self, ensembl_ids: list[str]) -> pl.DataFrame`. The returned frame has columns `ensembl_id`, `uniprot_id`, `primary_gene_name`, `ncbi_description`, `uniprot_description`, `subcellular_localisation`, all String; list-valued gget fields are joined with "; ".
  - Production defaults:
    - `ensembl_rest_symbol_lookup`: GET `https://rest.ensembl.org/lookup/symbol/homo_sapiens/{symbol}?content-type=application/json`, via `call_with_retries` with `retry_on=GGET_TRANSIENT_ERRORS`. A 400 response means unresolved and returns None.
    - `gget_info_lookup`: wraps `fetch_gget_info` and converts the result to polars.
  - `ColocProteinInfoTask` fields: `meta`, `display_task`, `max_loci: int`, `symbol_lookup: SymbolLookup`, `info_lookup: InfoLookup`.
  - Output parquet columns:
    - from the display table: `gene`, `protein`, `locus`, `pp_h4`;
    - `gene_symbol`, one per split gene;
    - then the InfoLookup columns.

- [ ] **Step 1: Write the failing tests.** They use fake lookups, which are plain classes that record their calls.
  - `test_one_row_per_gene_of_top_loci`: the display has 2 loci, `max_loci=1`, and the fake info returns a row per id. The output covers only the top locus's proteins.
  - `test_multi_gene_complex_split`: gene "IL12A_IL12B" gives 2 rows, with gene_symbol "IL12A" and "IL12B" and the same protein.
  - `test_unresolved_symbol_keeps_row_with_null_annotations`: the fake symbol lookup returns None for "FOO". Its row exists with `ensembl_id` null, and the info lookup was not called with None.
  - `test_empty_display_gives_empty_table`: 0 rows, and the schema is still present.

- [ ] **Step 2: Run them and verify they fail.**

- [ ] **Step 3: Implement.**
  - `execute`:
    1. read the display;
    2. take the top loci rows;
    3. explode the split symbols;
    4. map the symbols;
    5. call the info lookup once on the distinct non-null ids;
    6. left-join on ensembl_id;
    7. write zstd parquet via `write_parquet_table`.
  - The meta is a `ResultTableMeta` like the display task's, with sub_dir "analysis".

- [ ] **Step 4: Run the tests and verify they pass.** Then run green (tee to a log), which also covers the FetchGGetInfoTask refactor.

- [ ] **Step 5: Commit:** `"Add coloc protein annotation task (Ensembl symbol lookup + gget info)"`.

### Task 6: Generator helpers, DecodeME instances, and run

**Files:**
- Modify: `mecfs_bio/asset_generator/protein_coloc_asset_generator.py`. Add:
  - `generate_coloc_locus_stack_plot_task(asset_id: str, group: ColocTaskGroup, max_loci: int = 10) -> ColocLocusStackPlotTask`, with `gene_info_task=MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED` and `genetic_map_task=GENETIC_MAP_HG38`;
  - `generate_coloc_protein_info_task(asset_id: str, group: ColocTaskGroup, max_loci: int = 10) -> ColocProteinInfoTask`.
- Modify: `mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/decode_me_ppp_coloc.py`. Add `DECODE_ME_PPP_COLOC_EXPLORATORY_2_STACK_PLOTS` (asset id `decode_me_gwas_1_exploratory_2_ppp_coloc_stack_plots`) and `DECODE_ME_PPP_COLOC_EXPLORATORY_2_PROTEIN_INFO` (asset id `decode_me_gwas_1_exploratory_2_ppp_coloc_protein_info`).
- Create: `experiments/claude/protein_coloc/run_decode_me_exploratory_2_followups.py`. It force-rebuilds the display task, then builds the two new tasks, using the repo's existing runner pattern (copy from `experiments/claude/protein_coloc/rerun_coloc_gnomad_v4_1.py`). Its tests are the build itself.

- [ ] **Step 1: Add the generators and instances.** Run green (tee to a log). This is wiring only, so there are no structural tests.

- [ ] **Step 2: Run.**
  - Command: `pixi r python experiments/claude/protein_coloc/run_decode_me_exploratory_2_followups.py 2>&1 | tee experiments/claude/logs/decode_me_exploratory_2_followups.log`
  - Expected:
    - the coloc task is NOT re-executed (check the log);
    - the display table has 8 rows with the new columns;
    - 2 locus figures exist (chr1 and chr12);
    - the protein-info table has 8 rows.

- [ ] **Step 3: Inspect.**
  - Open both PNGs and confirm, for each locus:
    - the trait peak, the hit-protein peaks and the dotted trait-lead line line up;
    - the significance line is on the trait row only;
    - the recombination track is present.
  - Print the display table and the info table, and report to the user:
    - `same_sign` for each hit;
    - `n_snp`;
    - the UniProt one-liners.

- [ ] **Step 4: Commit** the generator, the DecodeME asset file and the run script: `"Add DecodeME exploratory PPP coloc stack plots and protein info"`.
