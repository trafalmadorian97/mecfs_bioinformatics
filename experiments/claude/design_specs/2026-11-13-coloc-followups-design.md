# Protein coloc follow-ups: R cross-check, richer display table, locus stack plots, protein info

Date: 2026-11-13. Branch: decode_me_coloc.

## Intent

The exploratory DecodeME x UKB-PPP coloc run (`DECODE_ME_PPP_COLOC_EXPLORATORY_2`: lead p < 1e-5,
max_pairs 2000, max_per_locus 30, MHC included) gave the first interesting ME/CFS results in a
while: 8 pairs with PP.H4 > 0.8 at two loci, all trans:

- chr1:7.94-8.94 Mb: DRAXIN, TCL1A, FCRL1.
- chr12:123.4-124.4 Mb: FGFBP1, CRIM1, FAM13A, LYPD3, NOTCH3.

The user wants to make sure these are correct, and to make them easy to inspect. In order of
priority:

0. Verify that our Python coloc.abf reproduces the reference R implementation on this run.
1. Add context columns to the display table.
2. Add stacked Manhattan plots per locus, to eyeball whether the trait and protein peaks coincide.
3. Add database information about the colocalizing proteins.

## Goals

- Show agreement (or a quantified disagreement) between our ColocAbf and R coloc::coloc.abf on the
  same inputs.
- Display table gains:
  - the trait lead p and the protein lead p;
  - whether the trait and protein effects have the same sign;
  - the number of variants passed to coloc.abf.
- One stacked figure per display-table locus, with these rows:
  - the trait Manhattan, with a genome-wide-significance line and the recombination rate;
  - one Manhattan per colocalizing protein;
  - gene tracks.
- A table of database annotations for the proteins at the top loci.
- Instances of all of these for `DECODE_ME_PPP_COLOC_EXPLORATORY_2`, built and inspected.

## Non-goals

- LD colouring from a reference LD matrix. This was dropped by the user: matching the GRCh37 Broad
  UKBB panel to GRCh38 windows has subtle pitfalls, and the download is large. Recombination rate
  stands in as rough LD context.
- Re-running TraitProteinColocTask. All new outputs derive from its existing result directory
  (raw_results.parquet and windows/*.parquet). No PPP or CSF rebuilds.
- Adding R coloc as a permanent dependency.

## Anti-goals

- Silently trusting a partial R comparison. If any pair disagrees beyond tolerance, that must be
  reported plainly, not averaged away.

## Constraints

- All commands via pixi. Tee `invoke green` output to a logfile. Never commit logs or PDFs. Never
  edit docs/.
- Display table column names stay short, so the table fits a docs page.
- Code changes do not invalidate the build cache, so the display table must be force-rebuilt.

## Approaches considered

- **LD colouring by the Broad UKBB LD matrix**, matched through the PPP variant index's hg19
  positions and downloaded and discarded on the fly. Rejected by the user (see Non-goals).
- **Recombination rate on the trait panel.** Chosen. The Eagle genetic map is already used for
  hg19 (GENETIC_MAP_HG19, plus the recombination track in polyfun_explain_plot_task.py). Eagle
  publishes the same table for hg38.
- **Extend ColocLocusPlotTask (one figure per hand-picked pair) versus a new per-locus task.** The
  new task chooses loci from the display table and stacks all colocalizing proteins in one figure.
  Shared helpers move to a common module. The existing task stays.
- **Sign recorded in the coloc task versus derived downstream.** Derived downstream, because
  recording it in the coloc task would force a coloc rerun. The saved windows already hold
  BETA_TRAIT and BETA_PROTEIN against the same effect allele.
- **R cross-check as a unit test versus an experiments script.** An experiments script, per the
  user. R coloc runs in a temporary environment
  (`pixi exec --spec r-coloc=5.2.3 --spec r-base ... Rscript`), so nothing is added to pyproject.

## Design

### 0. R cross-check script: experiments/claude/protein_coloc/r_coloc_crosscheck.py

Inputs:
- The exploratory_2 coloc result directory, located via the task's meta.
- The display table, for the list of hit pairs.

Steps:
- For each pair, export the saved window (BETA_TRAIT, SE_TRAIT, BETA_PROTEIN, SE_PROTEIN, EAF, N)
  to CSV in the scratch or experiment dir.
- An R script calls `coloc.abf`. Dataset 1 is the trait:

      type="cc", beta, varbeta = se^2, s = case fraction, snp

  Dataset 2 is the protein:

      type="quant", beta, varbeta, sdY = 1

  The priors are p1, p2 and p12 from the raw results row. The script writes PP.H0..H4.
- Python compares the R output to the raw_results PP columns.

Scope and tolerance:
- Run on every tested pair of the run (517 pairs; it is cheap), and report the display-table hits
  separately.
- Tolerance: max |dPP| < 1e-6. Report the max |dPP| overall and per hypothesis.
- The script must mirror whatever trait and protein models the task actually used:
  - CaseControlTraitModel(case_fraction) maps to cc, W = 0.2.
  - PROTEIN_TRAIT_MODEL = QuantitativeTraitModel(sd_y=1.0) maps to quant, sdY = 1.
- The builder checks that R's defaults do not differ for these settings: whether R needs MAF for
  type cc, and whether it uses s for cc when beta and varbeta are given.

Output and caveats:
- Output is a log in experiments/claude/logs (uncommitted) and a short summary appended to
  experiments/claude/protein_coloc/.
- This verifies the Bayes-factor and posterior arithmetic only. It does not verify
  harmonization or alignment upstream of the windows. The sign column and the plots help there.

### 1. Display table columns

New columns appended to DISPLAY_SCHEMA:

| Name | Type | Source |
|---|---|---|
| trait_p | Float64 | 10^-LEAD_MLOG10P |
| prot_p | Float64 | 10^-PROTEIN_LEAD_MLOG10P |
| same_sign | Boolean | sign(BETA_TRAIT) == sign(BETA_PROTEIN) at the protein's lead variant |
| n_snp | Int32 | N_ALIGNED |

- trait_p is the locus lead p over the whole trait window.
- prot_p is the protein lead p within the aligned window.
- For same_sign, the protein's lead variant is the max |BETA_PROTEIN/SE_PROTEIN| in the pair's
  saved window. The window is read from coloc_dir/windows/<pair_file_stem>.parquet, only for
  rows that pass the h4 filter.
- Keep the floats as floats. Docs formatting is the docs macro's job, so no string formatting
  goes in the parquet.
- Update the module docstring.
- Add a unit test: the existing fixture gets a windows dir, and the test checks same_sign for
  a positive/positive pair and a positive/negative pair.

### 2. ColocLocusStackPlotTask: mecfs_bio/build_system/task/protein_coloc/

Fields:
- coloc_task: TraitProteinColocTask. This gives the trait task and the result dir.
- display_task: ColocDisplayTableTask.
- gene_info_task: the GRCh38 gene table, the same one ColocLocusPlotTask uses.
- genetic_map_task: GENETIC_MAP_HG38.
- max_loci: int, default 10.

Loci and output:
- The loci are the distinct display-table loci, ranked by best PP.H4, keeping the first max_loci.
- The output is a directory with one PNG and one SVG per locus, plus plot_index.parquet
  (LOCUS_ID, N_PROTEINS, PNG, SVG).
- An empty display table is valid: it produces an empty index and no figures.

Rows, sharing a GRCh38 x-axis over the locus window:
1. Trait. -log10 p computed from beta/se of the full harmonized trait GWAS
   (coloc_task.trait_task) within the window, not the per-pair aligned subset. The panel has:
   - a dashed genome-wide-significance line (GWAS_SIGNIFICANCE_MLOG10P);
   - the lead variant marked;
   - recombination rate (cM/Mb) on a twin y-axis.
2. One row per display-table protein at the locus, in PP.H4 order, drawn from its saved window.
   Each row has:
   - the lead variant marked;
   - a vertical dotted line at the trait lead position, so offsets are visible;
   - a panel label with the gene, cis/trans status and PP.H4.
3. Gene tracks, using plot_gene_tracks.

New asset GENETIC_MAP_HG38:
- It sits beside GENETIC_MAP_HG19 and mirrors it, using the Eagle `genetic_map_hg38_withX.txt.gz`.
- Its md5 is recorded after the first download.
- Chromosome X is coded as 23 in Eagle; this does not matter for the autosomal loci here.

Code reuse:
- Hoist `_load_recomb` and `_plot_recomb` out of polyfun_explain_plot_task.py into a shared
  plotting module, with the docstring generalized beyond hg19.
- Reuse `_plot_association` and `_plot_genes` from coloc_locus_plot_task.py by moving them to
  that module or making them public.

Tests: a Task-level unit test with a synthetic coloc dir, display table, genes and map, checking:
- the figure files and index rows exist per locus;
- max_loci truncates.

### 3. ColocProteinInfoTask

Fields: display_task, max_loci (default 10, same ranking as above).

Steps:
- Take the distinct gene symbols of the display rows at the top loci.
- Split multi-gene complexes on "_", as in IL12A_IL12B.
- Map symbols to Ensembl gene IDs via the Ensembl REST symbol lookup, through call_with_retries
  (`GET /lookup/symbol/homo_sapiens/{symbol}`; or gget.search, at the builder's call).
- Call gget.info through call_with_retries, reusing what FetchGGetInfoTask does; refactor to
  share a helper rather than copy.

Output: a parquet with one row per (protein, gene) and these columns:
- the gene and protein display columns;
- locus;
- pp_h4;
- ensembl_id;
- uniprot_id;
- primary_gene_name;
- ncbi_description;
- uniprot_description;
- subcellular_localisation.

Symbols that fail to map get null annotation columns, with a warning log rather than an
abort.

Tests:
- A Task-level test with injected fake lookup and info callables (no network). This means the
  task takes them as fields with production defaults, per the no-monkeypatch convention.

### 4. Instances and run

Instances in mecfs_bio/assets/gwas/me_cfs/decode_me/analysis/decode_me_ppp_coloc.py, for
DECODE_ME_PPP_COLOC_EXPLORATORY_2:
- the plot task;
- the info task.

The display table instance already exists via the generator. Add a generator helper
next to generate_coloc_locus_plot_task, so other traits can reuse it.

Running:
- Force-rebuild the display table, then build the plot and info tasks.
- Inspect the figures for the chr1 and chr12 loci.

## Left to the builder

- Exact figure sizing and fonts. Follow ColocLocusPlotTask.
- Ensembl REST versus gget.search for the symbol mapping.
- Where the shared plotting helpers live.
