# Polyfun Annotation Weights Layer (Spec 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Produce a durable `gamma_c` annotation-weights asset by ridge-regressing the polyfun precomputed prior (`snpvar_bin`) on the 187 baseline-LF 2.2.UKB annotations.

**Architecture:** Three new Tasks + curated metadata. (1) reuse `DownloadFileTask` to fetch the 11GB baseline-LF tarball; (2) `BuildBaselineLFAnnotationParquetTask` extracts the per-chromosome `*.annot.gz` members and streams them into one (CHR,BP)-sorted annotation parquet; (3) `RidgeAnnotationWeightsTask` streams that parquet one chromosome at a time, accumulates per-chromosome cross-product sufficient statistics, picks alpha by leave-one-chromosome-out, and writes a weights directory (`weights.parquet` + `diagnostics.json`). A rule-based annotation→family classifier supports the hybrid attribution used in Spec 2.

**Tech Stack:** Python 3.13, polars, numpy, pyarrow, the repo build-system `Task` abstraction, pytest. All commands via `pixi r`.

**Spec:** `experiments/claude/design_specs/2026-08-25-polyfun-annotation-weights-design.md`

## Global Constraints

- Run everything via `pixi r <command>` (e.g. `pixi r pytest ...`, `pixi r python ...`).
- After significant changes run `pixi r invoke green` (capture to a logfile; `--testmon` may report "no tests ran").
- Tests: dependency injection only — NO monkeypatch/mock. Inject fake `WF`/`fetch`/dep-tasks.
- Prefer polars over pandas for new dataframe code.
- Path objects for filesystem paths; `PurePath` for relative paths.
- Use the repo subprocess wrapper `execute_command` (`mecfs_bio/util/subproc/run_command.py`), not raw subprocess.
- Reference column names via module-level constants, not repeated string literals.
- Write Task-level unit tests (the Task is the build system's public API); do not assert on log/error-message text.
- A Task's `create()` DERIVES its output meta from its dependency's meta (reuse `group`/`sub_group`/`sub_folder` from the source `ReferenceFileMeta`; `raise` on an unknown meta) — the `CompressedCSVToParquetTask.create` pattern. Never hardcode meta fields that the dependency already carries.
- Every function (implementation and test helpers) has BOTH parameter and return type annotations.
- Docstrings: no backticks around inline code, no RST.
- The regression target `snpvar_bin` is the polyfun *binned* prior actually used, not the original S-LDSC tau_c — reflect this in docstrings, do not call the coefficients tau_c.

## File Structure

New task package `mecfs_bio/build_system/task/annotation_weights/`:
- `__init__.py`
- `build_baseline_lf_annotation_parquet_task.py` — `BuildBaselineLFAnnotationParquetTask` + extraction/streaming helpers + `ANNOT_KEY_COLUMNS` constant.
- `ridge_annotation_weights_task.py` — `RidgeAnnotationWeightsTask` + Gram/LOCO helpers + weights column-name constants + `WEIGHTS_PARQUET_FILENAME`, `DIAGNOSTICS_JSON_FILENAME`.

New assets package `mecfs_bio/assets/reference_data/polyfun/annotations/`:
- `__init__.py`
- `baseline_lf_annotation_names.py` — `BASELINE_LF_ANNOTATION_NAMES: list[str]` (the 187 names).
- `annotation_families.py` — `AnnotationFamily` Literal + `family_for_annotation(name)`.
- `baseline_lf_annotations.py` — wiring: tarball `DownloadFileTask` instance + `BuildBaselineLFAnnotationParquetTask` instance.
- `annotation_ridge_weights.py` — wiring: `RidgeAnnotationWeightsTask` instance.

New tests under `test_mecfs_bio/unit/`:
- `assets/reference_data/polyfun/annotations/test_annotation_families.py`
- `build_system/task/annotation_weights/test_build_baseline_lf_annotation_parquet_task.py`
- `build_system/task/annotation_weights/test_ridge_annotation_weights_task.py`

(Create `__init__.py` in any new test directories that need them, matching sibling test dirs.)

---

### Task 1: Annotation names constant + family classifier

**Files:**
- Create: `mecfs_bio/assets/reference_data/polyfun/annotations/__init__.py` (empty)
- Create: `mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_annotation_names.py`
- Create: `mecfs_bio/assets/reference_data/polyfun/annotations/annotation_families.py`
- Test: `test_mecfs_bio/unit/assets/reference_data/polyfun/annotations/test_annotation_families.py` (+ `__init__.py` files up the new test path as needed)

**Interfaces:**
- Produces: `BASELINE_LF_ANNOTATION_NAMES: list[str]` (length 187, unique, order = annotation-file column order).
- Produces: `AnnotationFamily = Literal["non_synonymous","coding","conserved","promoter_or_enhancer","histone_marks","repressed","open_chromatin","maf_bins","ld_related_continuous","molecular_qtl","other"]` (11 families; taxonomy grounded in Weissbrod 2020's Supplementary Note grouping + Gazal 2017 / Hormozdiari 2018 for the groups it omits; `open_chromatin` is our one documented deviation from polyfun's "others" lumping).
- Produces: `family_for_annotation(name: str) -> AnnotationFamily`

- [ ] **Step 1: Create the names constant file**

Create `baseline_lf_annotation_names.py` with exactly this content (the 187 baseline-LF 2.2.UKB annotation column names, in file-column order):

```python
"""The 187 baseline-LF 2.2.UKB annotation column names, in annotation-file order.

Generated once from the header of a baselineLF2.2.UKB.<chr>.annot.gz file
(columns other than CHR, BP, SNP, CM). Committed as a constant so the family
map and its test do not depend on the 11GB annotation download.
"""

BASELINE_LF_ANNOTATION_NAMES: list[str] = [
    "Coding_UCSC_lowfreq",
    "Coding_UCSC_common",
    "Coding_UCSC.flanking.500_lowfreq",
    "Coding_UCSC.flanking.500_common",
    "Conserved_LindbladToh_lowfreq",
    "Conserved_LindbladToh_common",
    "Conserved_LindbladToh.flanking.500_lowfreq",
    "Conserved_LindbladToh.flanking.500_common",
    "CTCF_Hoffman_lowfreq",
    "CTCF_Hoffman_common",
    "CTCF_Hoffman.flanking.500_lowfreq",
    "CTCF_Hoffman.flanking.500_common",
    "DGF_ENCODE_lowfreq",
    "DGF_ENCODE_common",
    "DGF_ENCODE.flanking.500_lowfreq",
    "DGF_ENCODE.flanking.500_common",
    "DHS_peaks_Trynka_lowfreq",
    "DHS_peaks_Trynka_common",
    "DHS_Trynka_lowfreq",
    "DHS_Trynka_common",
    "DHS_Trynka.flanking.500_lowfreq",
    "DHS_Trynka.flanking.500_common",
    "Enhancer_Andersson_lowfreq",
    "Enhancer_Andersson_common",
    "Enhancer_Andersson.flanking.500_lowfreq",
    "Enhancer_Andersson.flanking.500_common",
    "Enhancer_Hoffman_lowfreq",
    "Enhancer_Hoffman_common",
    "Enhancer_Hoffman.flanking.500_lowfreq",
    "Enhancer_Hoffman.flanking.500_common",
    "FetalDHS_Trynka_lowfreq",
    "FetalDHS_Trynka_common",
    "FetalDHS_Trynka.flanking.500_lowfreq",
    "FetalDHS_Trynka.flanking.500_common",
    "H3K27ac_Hnisz_lowfreq",
    "H3K27ac_Hnisz_common",
    "H3K27ac_Hnisz.flanking.500_lowfreq",
    "H3K27ac_Hnisz.flanking.500_common",
    "H3K27ac_PGC2_lowfreq",
    "H3K27ac_PGC2_common",
    "H3K27ac_PGC2.flanking.500_lowfreq",
    "H3K27ac_PGC2.flanking.500_common",
    "H3K4me1_peaks_Trynka_lowfreq",
    "H3K4me1_peaks_Trynka_common",
    "H3K4me1_Trynka_lowfreq",
    "H3K4me1_Trynka_common",
    "H3K4me1_Trynka.flanking.500_lowfreq",
    "H3K4me1_Trynka.flanking.500_common",
    "H3K4me3_peaks_Trynka_lowfreq",
    "H3K4me3_peaks_Trynka_common",
    "H3K4me3_Trynka_lowfreq",
    "H3K4me3_Trynka_common",
    "H3K4me3_Trynka.flanking.500_lowfreq",
    "H3K4me3_Trynka.flanking.500_common",
    "H3K9ac_peaks_Trynka_lowfreq",
    "H3K9ac_peaks_Trynka_common",
    "H3K9ac_Trynka_lowfreq",
    "H3K9ac_Trynka_common",
    "H3K9ac_Trynka.flanking.500_lowfreq",
    "H3K9ac_Trynka.flanking.500_common",
    "Intron_UCSC_lowfreq",
    "Intron_UCSC_common",
    "Intron_UCSC.flanking.500_lowfreq",
    "Intron_UCSC.flanking.500_common",
    "PromoterFlanking_Hoffman_lowfreq",
    "PromoterFlanking_Hoffman_common",
    "PromoterFlanking_Hoffman.flanking.500_lowfreq",
    "PromoterFlanking_Hoffman.flanking.500_common",
    "Promoter_UCSC_lowfreq",
    "Promoter_UCSC_common",
    "Promoter_UCSC.flanking.500_lowfreq",
    "Promoter_UCSC.flanking.500_common",
    "Repressed_Hoffman_lowfreq",
    "Repressed_Hoffman_common",
    "Repressed_Hoffman.flanking.500_lowfreq",
    "Repressed_Hoffman.flanking.500_common",
    "SuperEnhancer_Hnisz_lowfreq",
    "SuperEnhancer_Hnisz_common",
    "SuperEnhancer_Hnisz.flanking.500_lowfreq",
    "SuperEnhancer_Hnisz.flanking.500_common",
    "TFBS_ENCODE_lowfreq",
    "TFBS_ENCODE_common",
    "TFBS_ENCODE.flanking.500_lowfreq",
    "TFBS_ENCODE.flanking.500_common",
    "Transcr_Hoffman_lowfreq",
    "Transcr_Hoffman_common",
    "Transcr_Hoffman.flanking.500_lowfreq",
    "Transcr_Hoffman.flanking.500_common",
    "TSS_Hoffman_lowfreq",
    "TSS_Hoffman_common",
    "TSS_Hoffman.flanking.500_lowfreq",
    "TSS_Hoffman.flanking.500_common",
    "UTR_3_UCSC_lowfreq",
    "UTR_3_UCSC_common",
    "UTR_3_UCSC.flanking.500_lowfreq",
    "UTR_3_UCSC.flanking.500_common",
    "UTR_5_UCSC_lowfreq",
    "UTR_5_UCSC_common",
    "UTR_5_UCSC.flanking.500_lowfreq",
    "UTR_5_UCSC.flanking.500_common",
    "WeakEnhancer_Hoffman_lowfreq",
    "WeakEnhancer_Hoffman_common",
    "WeakEnhancer_Hoffman.flanking.500_lowfreq",
    "WeakEnhancer_Hoffman.flanking.500_common",
    "GERP.NS_lowfreq",
    "GERP.NS_common",
    "GERP.RSsup4_lowfreq",
    "GERP.RSsup4_common",
    "MAFbin_lowfreq_1",
    "MAFbin_lowfreq_2",
    "MAFbin_lowfreq_3",
    "MAFbin_lowfreq_4",
    "MAFbin_lowfreq_5",
    "MAFbin_lowfreq_6",
    "MAFbin_lowfreq_7",
    "MAFbin_lowfreq_8",
    "MAFbin_lowfreq_9",
    "MAFbin_lowfreq_10",
    "MAFbin_frequent_1",
    "MAFbin_frequent_2",
    "MAFbin_frequent_3",
    "MAFbin_frequent_4",
    "MAFbin_frequent_5",
    "MAFbin_frequent_6",
    "MAFbin_frequent_7",
    "MAFbin_frequent_8",
    "MAFbin_frequent_9",
    "MAFbin_frequent_10",
    "MAF_Adj_Predicted_Allele_Age_common",
    "MAF_Adj_LLD_AFR_lowfreq",
    "MAF_Adj_LLD_AFR_common",
    "Recomb_Rate_10kb_lowfreq",
    "Recomb_Rate_10kb_common",
    "Nucleotide_Diversity_10kb_lowfreq",
    "Nucleotide_Diversity_10kb_common",
    "Backgrd_Selection_Stat_lowfreq",
    "Backgrd_Selection_Stat_common",
    "CpG_Content_50kb_lowfreq",
    "CpG_Content_50kb_common",
    "MAF_Adj_ASMC_lowfreq",
    "MAF_Adj_ASMC_common",
    "GTEx_eQTL_MaxCPP_common",
    "BLUEPRINT_H3K27acQTL_MaxCPP_common",
    "BLUEPRINT_H3K4me1QTL_MaxCPP_common",
    "BLUEPRINT_DNA_methylation_MaxCPP_common",
    "synonymous_lowfreq",
    "synonymous_common",
    "non_synonymous_lowfreq",
    "non_synonymous_common",
    "Conserved_Vertebrate_phastCons46way_lowfreq",
    "Conserved_Vertebrate_phastCons46way_common",
    "Conserved_Vertebrate_phastCons46way.flanking.500_lowfreq",
    "Conserved_Vertebrate_phastCons46way.flanking.500_common",
    "Conserved_Mammal_phastCons46way_lowfreq",
    "Conserved_Mammal_phastCons46way_common",
    "Conserved_Mammal_phastCons46way.flanking.500_lowfreq",
    "Conserved_Mammal_phastCons46way.flanking.500_common",
    "Conserved_Primate_phastCons46way_lowfreq",
    "Conserved_Primate_phastCons46way_common",
    "Conserved_Primate_phastCons46way.flanking.500_lowfreq",
    "Conserved_Primate_phastCons46way.flanking.500_common",
    "BivFlnk_lowfreq",
    "BivFlnk_common",
    "BivFlnk.flanking.500_lowfreq",
    "BivFlnk.flanking.500_common",
    "Human_Promoter_Villar_lowfreq",
    "Human_Promoter_Villar_common",
    "Human_Promoter_Villar.flanking.500_lowfreq",
    "Human_Promoter_Villar.flanking.500_common",
    "Human_Enhancer_Villar_lowfreq",
    "Human_Enhancer_Villar_common",
    "Human_Enhancer_Villar.flanking.500_lowfreq",
    "Human_Enhancer_Villar.flanking.500_common",
    "Ancient_Sequence_Age_Human_Promoter_lowfreq",
    "Ancient_Sequence_Age_Human_Promoter_common",
    "Ancient_Sequence_Age_Human_Promoter.flanking.500_lowfreq",
    "Ancient_Sequence_Age_Human_Promoter.flanking.500_common",
    "Ancient_Sequence_Age_Human_Enhancer_lowfreq",
    "Ancient_Sequence_Age_Human_Enhancer_common",
    "Ancient_Sequence_Age_Human_Enhancer.flanking.500_lowfreq",
    "Ancient_Sequence_Age_Human_Enhancer.flanking.500_common",
    "Human_Enhancer_Villar_Species_Enhancer_Count_lowfreq",
    "Human_Enhancer_Villar_Species_Enhancer_Count_common",
    "Human_Promoter_Villar_ExAC_lowfreq",
    "Human_Promoter_Villar_ExAC_common",
    "Human_Promoter_Villar_ExAC.flanking.500_lowfreq",
    "Human_Promoter_Villar_ExAC.flanking.500_common",
]
```

- [ ] **Step 2: Write the failing family-classifier test**

Create the test (and `__init__.py` files along the new test path). It covers every annotation, so it fails until the classifier exists:

```python
from typing import get_args

from mecfs_bio.assets.reference_data.polyfun.annotations.annotation_families import (
    AnnotationFamily,
    family_for_annotation,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotation_names import (
    BASELINE_LF_ANNOTATION_NAMES,
)


def test_names_constant_is_187_unique():
    assert len(BASELINE_LF_ANNOTATION_NAMES) == 187
    assert len(set(BASELINE_LF_ANNOTATION_NAMES)) == 187


def test_every_annotation_maps_to_a_valid_family():
    valid = set(get_args(AnnotationFamily))
    for name in BASELINE_LF_ANNOTATION_NAMES:
        assert family_for_annotation(name) in valid


def test_representative_family_assignments():
    cases = {
        "non_synonymous_lowfreq": "non_synonymous",
        "Coding_UCSC_common": "coding",
        "synonymous_common": "coding",
        "UTR_3_UCSC_common": "coding",
        "Conserved_Primate_phastCons46way_common": "conserved",
        "GERP.RSsup4_common": "conserved",
        "GERP.NS_common": "conserved",
        "Promoter_UCSC_common": "promoter_or_enhancer",
        "TSS_Hoffman_common": "promoter_or_enhancer",
        "SuperEnhancer_Hnisz_common": "promoter_or_enhancer",
        "BivFlnk.flanking.500_common": "promoter_or_enhancer",
        "H3K27ac_Hnisz_common": "histone_marks",
        "H3K4me1_peaks_Trynka_common": "histone_marks",
        "Repressed_Hoffman_common": "repressed",
        "DHS_Trynka_common": "open_chromatin",
        "DHS_peaks_Trynka_common": "open_chromatin",
        "FetalDHS_Trynka_common": "open_chromatin",
        "DGF_ENCODE_common": "open_chromatin",
        "TFBS_ENCODE_common": "other",
        "CTCF_Hoffman_common": "other",
        "Transcr_Hoffman_common": "other",
        "Intron_UCSC_common": "other",
        "MAFbin_frequent_3": "maf_bins",
        "CpG_Content_50kb_common": "ld_related_continuous",
        "Recomb_Rate_10kb_common": "ld_related_continuous",
        "Backgrd_Selection_Stat_common": "ld_related_continuous",
        "MAF_Adj_ASMC_common": "ld_related_continuous",
        "GTEx_eQTL_MaxCPP_common": "molecular_qtl",
        "BLUEPRINT_H3K27acQTL_MaxCPP_common": "molecular_qtl",
    }
    for name, family in cases.items():
        assert family_for_annotation(name) == family


def test_open_chromatin_membership_is_exactly_dhs_and_dgf():
    open_chrom = {
        n for n in BASELINE_LF_ANNOTATION_NAMES
        if family_for_annotation(n) == "open_chromatin"
    }
    # every open_chromatin member is a DHS/FetalDHS/DGF accessibility annotation
    assert all(("DHS" in n) or ("DGF" in n) for n in open_chrom)
    # and every DHS/DGF accessibility annotation is captured
    assert all(
        family_for_annotation(n) == "open_chromatin"
        for n in BASELINE_LF_ANNOTATION_NAMES
        if ("DHS" in n) or ("DGF" in n)
    )
```

- [ ] **Step 3: Run the test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/assets/reference_data/polyfun/annotations/test_annotation_families.py -v`
Expected: FAIL (import error — `annotation_families` not defined).

- [ ] **Step 4: Implement the classifier**

Create `annotation_families.py`:

```python
"""Map each baseline-LF 2.2.UKB annotation to one of eleven functional families.

Used for the hybrid attribution in the polyfun explainability pipeline: ridge
weights are fit on all 187 annotations, but contributions are aggregated to
families for headline reporting.

The family taxonomy is grounded in published sources, not invented:
  - The functional-group names (non_synonymous, coding, conserved,
    promoter_or_enhancer, histone_marks, repressed, other) are the grouping the
    polyfun authors themselves use for these annotations in the sub-additive
    simulation of their Supplementary Note (Weissbrod et al. 2020, Nat Genet).
  - maf_bins and ld_related_continuous are the MAF-bin and LD-related continuous
    annotation groups introduced in Gazal et al. 2017 (Nat Genet) baseline-LD
    (the Continuous rows of Gazal et al. 2018 Supplementary Table 1).
  - molecular_qtl are the MaxCPP fine-mapped molecular-QTL annotations of
    Hormozdiari et al. 2018 (Nat Genet).
  - open_chromatin (DHS/FetalDHS/DGF accessibility annotations) is the ONE
    deliberate refinement of polyfun's scheme, which otherwise lumps these into
    "others"; broken out so the explainability figure can show an accessibility
    panel. TFBS/CTCF/Transcribed/Intron remain in other, as in polyfun's scheme.

Per-annotation assignment is rule-based (keyword + a small override set) and
follows the annotation names and their source datasets (Gazal et al. 2018
Supplementary Table 1). The test in test_annotation_families.py asserts every one
of the 187 annotations resolves to a valid family.
"""

from typing import Literal

AnnotationFamily = Literal[
    "non_synonymous",
    "coding",
    "conserved",
    "promoter_or_enhancer",
    "histone_marks",
    "repressed",
    "open_chromatin",
    "maf_bins",
    "ld_related_continuous",
    "molecular_qtl",
    "other",
]

# Explicit overrides, matched as substrings and checked BEFORE the keyword rules.
# These are the continuous/special/molecular-QTL annotations whose family is not
# implied by a plain functional keyword (or that must beat a later keyword).
_OVERRIDES: tuple[tuple[str, AnnotationFamily], ...] = (
    # MaxCPP molecular-QTL (must win over the "H3K"/histone keyword) - Hormozdiari 2018
    ("MaxCPP", "molecular_qtl"),
    # LD-related continuous - Gazal 2017
    ("Backgrd_Selection", "ld_related_continuous"),
    ("Nucleotide_Diversity", "ld_related_continuous"),
    ("CpG_Content", "ld_related_continuous"),
    ("MAF_Adj_LLD_AFR", "ld_related_continuous"),
    ("Recomb_Rate", "ld_related_continuous"),
    ("MAF_Adj_ASMC", "ld_related_continuous"),
    ("Predicted_Allele_Age", "ld_related_continuous"),
    # GERP NS is a continuous conservation annotation -> conserved (by function)
    ("GERP.NS", "conserved"),
    # flanking bivalent TSS/enhancer -> promoter_or_enhancer
    ("BivFlnk", "promoter_or_enhancer"),
    # genic, non-accessibility -> other
    ("Intron_UCSC", "other"),
)

# Ordered keyword rules; the first family whose keyword appears in the name wins.
# "non_synonymous" is checked before "synonymous" so it is not swallowed by coding.
_KEYWORD_RULES: tuple[tuple[tuple[str, ...], AnnotationFamily], ...] = (
    (("MAFbin",), "maf_bins"),
    (("non_synonymous",), "non_synonymous"),
    (("synonymous", "Coding", "UTR_"), "coding"),
    (("Conserved", "phastCons", "GERP"), "conserved"),
    (("DHS", "DGF"), "open_chromatin"),
    (("Promoter", "TSS"), "promoter_or_enhancer"),
    (("Enhancer",), "promoter_or_enhancer"),
    (("H3K",), "histone_marks"),
    (("Repressed",), "repressed"),
    (("TFBS", "CTCF", "Transcr"), "other"),
)


def family_for_annotation(name: str) -> AnnotationFamily:
    """Return the functional family for a baseline-LF annotation column name."""
    for pattern, family in _OVERRIDES:
        if pattern in name:
            return family
    for keywords, family in _KEYWORD_RULES:
        if any(keyword in name for keyword in keywords):
            return family
    raise ValueError(f"No family rule matched annotation {name!r}")
```

- [ ] **Step 5: Run the test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/assets/reference_data/polyfun/annotations/test_annotation_families.py -v`
Expected: PASS (all three tests).

- [ ] **Step 6: Commit**

```bash
git add mecfs_bio/assets/reference_data/polyfun/annotations test_mecfs_bio/unit/assets/reference_data/polyfun/annotations
git commit -m "feat: baseline-LF annotation names + family classifier"
```

---

### Task 2: BuildBaselineLFAnnotationParquetTask

**Files:**
- Create: `mecfs_bio/build_system/task/annotation_weights/__init__.py` (empty)
- Create: `mecfs_bio/build_system/task/annotation_weights/build_baseline_lf_annotation_parquet_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_build_baseline_lf_annotation_parquet_task.py` (+ `__init__.py`)

**Interfaces:**
- Consumes: a tarball dependency Task producing a `FileAsset` pointing at a `baselineLF_v2.2.UKB.tar.gz`-shaped archive (members `<dir>/baselineLF2.2.UKB.<chr>.annot.gz`, LDSC-tab-separated: `CHR, BP, SNP, CM` + annotation columns; plus decoy non-annot members).
- Produces: `BuildBaselineLFAnnotationParquetTask` (frozen Task) with `.create(asset_id: str, tarball_task: Task) -> BuildBaselineLFAnnotationParquetTask`, emitting a single `FileAsset` parquet with columns `CHR (int), BP (int), SNP (str), CM (float), <annotations: f32...>`, deduped on `SNP`, ordered by `(CHR, BP)`.
- Produces: `ANNOT_KEY_COLUMNS: list[str] = ["CHR", "BP", "SNP", "CM"]`.

- [ ] **Step 1: Write the failing test**

```python
import gzip
import tarfile
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    BuildBaselineLFAnnotationParquetTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _write_annot_gz(path: Path, rows: list[dict]) -> None:
    header = ["CHR", "BP", "SNP", "CM", "annotA", "annotB"]
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(str(r[c]) for c in header))
    path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))


def _build_fake_tarball(tmp_path: Path) -> Path:
    src = tmp_path / "src" / "baselineLF_v2.2.UKB"
    src.mkdir(parents=True)
    # chr2 out of order relative to chr1 to prove global (CHR,BP) sort;
    # a duplicate rsid on chr1 to prove dedup.
    _write_annot_gz(
        src / "baselineLF2.2.UKB.1.annot.gz",
        [
            {"CHR": 1, "BP": 100, "SNP": "rs1", "CM": 0.1, "annotA": 1, "annotB": 0.5},
            {"CHR": 1, "BP": 200, "SNP": "rs2", "CM": 0.2, "annotA": 0, "annotB": 0.25},
            {"CHR": 1, "BP": 200, "SNP": "rs2", "CM": 0.2, "annotA": 0, "annotB": 0.25},
        ],
    )
    _write_annot_gz(
        src / "baselineLF2.2.UKB.2.annot.gz",
        [
            {"CHR": 2, "BP": 50, "SNP": "rs3", "CM": 0.3, "annotA": 1, "annotB": 0.75},
        ],
    )
    # decoy non-annotation member that must be ignored
    (src / "baselineLF2.2.UKB.1.l2.ldscore.gz").write_bytes(gzip.compress(b"junk\n"))
    tarball = tmp_path / "baselineLF_v2.2.UKB.tar.gz"
    with tarfile.open(tarball, "w:gz") as tar:
        tar.add(src, arcname="baselineLF_v2.2.UKB")
    return tarball


def test_builds_sorted_deduped_annotation_parquet(tmp_path: Path):
    tarball = _build_fake_tarball(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    tarball_task = FakeTask(
        SimpleFileMeta("annot_tarball", read_spec=None)
    )
    task = BuildBaselineLFAnnotationParquetTask.create(
        asset_id="annot_parquet", tarball_task=tarball_task
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "annot_tarball":
            return FileAsset(tarball)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    df = pl.read_parquet(result.path)

    # dedup: rs2 appears once -> 3 unique SNPs total
    assert df.height == 3
    assert df["SNP"].to_list() == ["rs1", "rs2", "rs3"]  # (CHR,BP) sorted
    assert df["CHR"].to_list() == [1, 1, 2]
    # annotation columns present and float32
    assert df.schema["annotA"] == pl.Float32
    assert df.schema["annotB"] == pl.Float32
    assert set(df.columns) == {"CHR", "BP", "SNP", "CM", "annotA", "annotB"}
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_build_baseline_lf_annotation_parquet_task.py -v`
Expected: FAIL (import error — task not defined).

- [ ] **Step 3: Implement the task**

Create `build_baseline_lf_annotation_parquet_task.py`:

```python
"""Build one sorted annotation parquet from the baseline-LF 2.2.UKB tarball.

Extracts the per-chromosome baselineLF2.2.UKB.<chr>.annot.gz members (LDSC text:
CHR, BP, SNP, CM + annotation columns), casts the annotation columns to float32,
dedups multiallelic rsid duplicates, and streams them, in chromosome order, into
a single (CHR, BP)-sorted parquet. The single sorted file gives the ridge weights
task a cheap scan and the downstream explainability tasks predicate-pushdown
per-locus lookups.
"""

import re
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    write_df_according_to_format,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

ANNOT_KEY_COLUMNS: list[str] = ["CHR", "BP", "SNP", "CM"]
_SNP_COL = "SNP"
_CHR_COL = "CHR"
_BP_COL = "BP"
_ANNOT_MEMBER_RE = re.compile(r"baselineLF2\.2\.UKB\.(\d+)\.annot\.gz$")


@frozen
class BuildBaselineLFAnnotationParquetTask(Task):
    meta: Meta
    tarball_task: Task

    @property
    def deps(self) -> list["Task"]:
        return [self.tarball_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        tarball_asset = fetch(self.tarball_task.asset_id)
        assert isinstance(tarball_asset, FileAsset)
        members = _list_annot_members(tarball_asset.path)
        extract_dir = scratch_dir / "annot_gz"
        extract_dir.mkdir(parents=True, exist_ok=True)
        _extract_members(tarball_asset.path, members, extract_dir)

        lazy_frames = []
        for chrom in sorted(members):
            member_path = extract_dir / members[chrom]
            lazy_frames.append(_scan_one_chromosome(member_path))
        combined = pl.concat(lazy_frames, how="vertical").sort([_CHR_COL, _BP_COL])

        out_path = scratch_dir / self.meta.asset_id
        write_df_according_to_format(
            df=combined, out_path=out_path, out_format=ParquetOutFormat()
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, asset_id: str, tarball_task: Task
    ) -> "BuildBaselineLFAnnotationParquetTask":
        # Derive the output meta from the tarball dependency's meta (reuse
        # group/sub_group/sub_folder), the CompressedCSVToParquetTask.create
        # pattern. Only the id, extension and read_spec differ.
        source_meta = tarball_task.meta
        if not isinstance(source_meta, ReferenceFileMeta):
            raise ValueError(f"Unknown meta for tarball task: {source_meta}")
        meta = ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
            extension=".parquet",
            read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
        )
        return cls(meta=meta, tarball_task=tarball_task)


def _list_annot_members(tarball: Path) -> dict[int, str]:
    """Map chromosome -> archive member path for every *.annot.gz member."""
    listing = execute_command(["tar", "-tzf", str(tarball)])
    members: dict[int, str] = {}
    for line in listing.splitlines():
        name = line.strip()
        match = _ANNOT_MEMBER_RE.search(name)
        if match:
            members[int(match.group(1))] = name
    assert members, f"no *.annot.gz members found in {tarball}"
    return members


def _extract_members(tarball: Path, members: dict[int, str], dest: Path) -> None:
    execute_command(
        ["tar", "-xzf", str(tarball), "-C", str(dest), *members.values()]
    )


def _scan_one_chromosome(member_path: Path) -> pl.LazyFrame:
    lazy = pl.scan_csv(
        member_path, separator="\t", infer_schema_length=None
    )
    schema = lazy.collect_schema()
    annot_cols = [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]
    return (
        lazy.with_columns([pl.col(c).cast(pl.Float32) for c in annot_cols])
        .unique(subset=_SNP_COL, keep="first", maintain_order=True)
    )
```

Note on the extracted member path: `tar -xzf ... -C dest <member>` writes to `dest/<member>` preserving the archive's internal directory (e.g. `dest/baselineLF_v2.2.UKB/baselineLF2.2.UKB.1.annot.gz`), which is exactly the member string returned by `_list_annot_members`, so `extract_dir / members[chrom]` resolves correctly.

- [ ] **Step 4: Run the test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_build_baseline_lf_annotation_parquet_task.py -v`
Expected: PASS. (If polars streaming over gzip errors, the fallback is `pl.read_csv(...).lazy()` per chromosome — still sunk through the same writer — but try the scan path first.)

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/build_system/task/annotation_weights test_mecfs_bio/unit/build_system/task/annotation_weights
git commit -m "feat: BuildBaselineLFAnnotationParquetTask"
```

---

### Task 3: RidgeAnnotationWeightsTask

**Files:**
- Create: `mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py`

**Interfaces:**
- Consumes: an annotation-parquet dep Task (`FileAsset` parquet: `CHR, BP, SNP, CM, <annotations>`) and a snpvar-meta dep Task (`FileAsset` parquet with columns `SNP`, `snpvar_bin`).
- Produces: `RidgeAnnotationWeightsTask` (frozen Task) with `.create(asset_id, annotation_parquet_task, snpvar_meta_task, alphas=(...)) -> RidgeAnnotationWeightsTask`, emitting a `DirectoryAsset` containing `weights.parquet` (columns `annotation, gamma_raw, gamma_standardized, family`) and `diagnostics.json` (`alpha, intercept, heldout_r2_per_chrom, mean_heldout_r2, n_variants`).
- Produces module constants: `WEIGHTS_PARQUET_FILENAME = "weights.parquet"`, `DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"`, `ANNOTATION_COL = "annotation"`, `GAMMA_RAW_COL = "gamma_raw"`, `GAMMA_STANDARDIZED_COL = "gamma_standardized"`, `FAMILY_COL = "family"`, `SNPVAR_COL = "snpvar_bin"`.

- [ ] **Step 1: Write the failing test**

```python
import json
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    DIAGNOSTICS_JSON_FILENAME,
    GAMMA_RAW_COL,
    WEIGHTS_PARQUET_FILENAME,
    RidgeAnnotationWeightsTask,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf


def _make_inputs(tmp_path: Path) -> tuple[Path, Path, dict[str, float]]:
    rng = np.random.default_rng(0)
    n_per_chrom = 400
    rows: list[dict] = []
    truth = {"annotA": 2.0, "annotB": -1.0, "annotC": 0.5}
    for chrom in (1, 2):
        a = rng.normal(size=n_per_chrom)
        b = rng.normal(size=n_per_chrom)
        c = rng.normal(size=n_per_chrom)
        y = 3.0 + truth["annotA"] * a + truth["annotB"] * b + truth["annotC"] * c
        for i in range(n_per_chrom):
            rows.append(
                {
                    "CHR": chrom,
                    "BP": i + 1,
                    "SNP": f"rs{chrom}_{i}",
                    "CM": 0.0,
                    "annotA": a[i],
                    "annotB": b[i],
                    "annotC": c[i],
                    "snpvar_bin": y[i],
                }
            )
    frame = pl.DataFrame(rows)
    annot_path = tmp_path / "annot.parquet"
    frame.drop("snpvar_bin").write_parquet(annot_path)
    meta_path = tmp_path / "meta.parquet"
    frame.select("SNP", "snpvar_bin").write_parquet(meta_path)
    return annot_path, meta_path, truth


def test_recovers_known_linear_weights(tmp_path: Path):
    annot_path, meta_path, truth = _make_inputs(tmp_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()

    annot_task = FakeTask(
        SimpleFileMeta("annot", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )
    meta_task = FakeTask(
        SimpleFileMeta("meta", read_spec=DataFrameReadSpec(DataFrameParquetFormat()))
    )
    task = RidgeAnnotationWeightsTask.create(
        asset_id="ridge_weights",
        annotation_parquet_task=annot_task,
        snpvar_meta_task=meta_task,
        alphas=(1e-4, 1e-2, 1.0, 100.0),
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "annot":
            return FileAsset(annot_path)
        if asset_id == "meta":
            return FileAsset(meta_path)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)

    weights = pl.read_parquet(result.path / WEIGHTS_PARQUET_FILENAME)
    got = dict(zip(weights[ANNOTATION_COL], weights[GAMMA_RAW_COL]))
    for name, expected in truth.items():
        assert abs(got[name] - expected) < 1e-2

    diagnostics = json.loads((result.path / DIAGNOSTICS_JSON_FILENAME).read_text())
    assert diagnostics["mean_heldout_r2"] > 0.999
    assert diagnostics["n_variants"] == 800
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py -v`
Expected: FAIL (import error).

- [ ] **Step 3: Implement the task**

Create `ridge_annotation_weights_task.py`:

```python
"""Fit a ridge surrogate of the polyfun prior on the baseline-LF annotations.

Regresses snpvar_bin (the polyfun binned per-SNP heritability prior actually used
in fine mapping, not the original S-LDSC tau_c) on the 187 annotations, genome
wide. The fit is done from per-chromosome cross-product sufficient statistics, so
the full design matrix is never held in memory; alpha is chosen by
leave-one-chromosome-out. Outputs raw-scale coefficients gamma_raw (used by the
explainability contrast) and standardized coefficients gamma_standardized (for
global importance ranking), plus each annotation's family.
"""

import json
from pathlib import Path

import numpy as np
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.assets.reference_data.polyfun.annotations.annotation_families import (
    family_for_annotation,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    ANNOT_KEY_COLUMNS,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

WEIGHTS_PARQUET_FILENAME = "weights.parquet"
DIAGNOSTICS_JSON_FILENAME = "diagnostics.json"
ANNOTATION_COL = "annotation"
GAMMA_RAW_COL = "gamma_raw"
GAMMA_STANDARDIZED_COL = "gamma_standardized"
FAMILY_COL = "family"
SNPVAR_COL = "snpvar_bin"
_CHR_COL = "CHR"
_SNP_COL = "SNP"

_DEFAULT_ALPHAS: tuple[float, ...] = (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


@frozen
class _ChromStats:
    """Cross-product sufficient statistics for one chromosome (standardized-free,
    raw annotation scale). p = number of annotations."""

    n: int
    sx: np.ndarray  # (p,) sum of annotations
    sxx: np.ndarray  # (p, p) sum of a a^T
    sxy: np.ndarray  # (p,) sum of a * y
    sy: float
    syy: float


@frozen
class RidgeAnnotationWeightsTask(Task):
    meta: Meta
    annotation_parquet_task: Task
    snpvar_meta_task: Task
    alphas: tuple[float, ...] = _DEFAULT_ALPHAS

    @property
    def deps(self) -> list["Task"]:
        return [self.annotation_parquet_task, self.snpvar_meta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        annot_asset = fetch(self.annotation_parquet_task.asset_id)
        assert isinstance(annot_asset, FileAsset)
        meta_asset = fetch(self.snpvar_meta_task.asset_id)
        meta = (
            scan_dataframe_asset(meta_asset, self.snpvar_meta_task.meta)
            .select(_SNP_COL, SNPVAR_COL)
            .collect()
            .to_polars()
            .unique(subset=_SNP_COL)
        )

        annot_columns = _annotation_columns(annot_asset.path)
        per_chrom = _accumulate_per_chromosome(
            annot_asset.path, annot_columns, meta
        )
        alpha, mean_r2, r2_per_chrom = _select_alpha_loco(per_chrom, self.alphas)
        gamma_raw, gamma_std, intercept = _fit_all(per_chrom, alpha)

        weights = pl.DataFrame(
            {
                ANNOTATION_COL: annot_columns,
                GAMMA_RAW_COL: gamma_raw,
                GAMMA_STANDARDIZED_COL: gamma_std,
                FAMILY_COL: [family_for_annotation(c) for c in annot_columns],
            }
        )
        weights.write_parquet(scratch_dir / WEIGHTS_PARQUET_FILENAME)
        diagnostics = {
            "alpha": alpha,
            "intercept": intercept,
            "heldout_r2_per_chrom": r2_per_chrom,
            "mean_heldout_r2": mean_r2,
            "n_variants": int(sum(s.n for s in per_chrom.values())),
        }
        (scratch_dir / DIAGNOSTICS_JSON_FILENAME).write_text(
            json.dumps(diagnostics, indent=2, sort_keys=True)
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        annotation_parquet_task: Task,
        snpvar_meta_task: Task,
        alphas: tuple[float, ...] = _DEFAULT_ALPHAS,
    ) -> "RidgeAnnotationWeightsTask":
        # Derive the output directory meta from the primary dependency (the
        # annotation parquet), reusing group/sub_group/sub_folder - the
        # CompressedCSVToParquetTask.create pattern applied to a DirMeta.
        source_meta = annotation_parquet_task.meta
        if not isinstance(source_meta, ReferenceFileMeta):
            raise ValueError(
                f"Unknown meta for annotation parquet task: {source_meta}"
            )
        meta = ReferenceDataDirectoryMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            id=AssetId(asset_id),
        )
        return cls(
            meta=meta,
            annotation_parquet_task=annotation_parquet_task,
            snpvar_meta_task=snpvar_meta_task,
            alphas=alphas,
        )


def _annotation_columns(annot_path: Path) -> list[str]:
    schema = pl.scan_parquet(annot_path).collect_schema()
    return [c for c in schema.names() if c not in ANNOT_KEY_COLUMNS]


def _accumulate_per_chromosome(
    annot_path: Path, annot_columns: list[str], meta: pl.DataFrame
) -> dict[int, _ChromStats]:
    chroms = (
        pl.scan_parquet(annot_path)
        .select(_CHR_COL)
        .unique()
        .collect()
        .to_series()
        .to_list()
    )
    per_chrom: dict[int, _ChromStats] = {}
    for chrom in sorted(chroms):
        frame = (
            pl.scan_parquet(annot_path)
            .filter(pl.col(_CHR_COL) == chrom)
            .collect()
            .to_polars()
            .join(meta, on=_SNP_COL, how="inner")
        )
        x = frame.select(annot_columns).to_numpy().astype(np.float64)
        y = frame.select(SNPVAR_COL).to_numpy().ravel().astype(np.float64)
        per_chrom[chrom] = _ChromStats(
            n=x.shape[0],
            sx=x.sum(0),
            sxx=x.T @ x,
            sxy=x.T @ y,
            sy=float(y.sum()),
            syy=float(y @ y),
        )
    return per_chrom


def _combine(stats: list[_ChromStats]) -> _ChromStats:
    """Sum per-chromosome sufficient statistics into one aggregate."""
    return _ChromStats(
        n=sum(s.n for s in stats),
        sx=sum((s.sx for s in stats), start=np.zeros_like(stats[0].sx)),
        sxx=sum((s.sxx for s in stats), start=np.zeros_like(stats[0].sxx)),
        sxy=sum((s.sxy for s in stats), start=np.zeros_like(stats[0].sxy)),
        sy=sum(s.sy for s in stats),
        syy=sum(s.syy for s in stats),
    )


@frozen
class _StandardizedSystem:
    """The centered+standardized ridge system for a set of annotations."""

    g_std: np.ndarray  # (p, p) standardized Gram
    b_std: np.ndarray  # (p,) standardized cross-term
    mean: np.ndarray  # (p,) per-annotation mean
    sd: np.ndarray  # (p,) per-annotation std (zeros replaced by 1)
    mean_y: float


def _standardized_system(stats: _ChromStats) -> _StandardizedSystem:
    """Build the centered+standardized ridge system from raw cross-products."""
    n = stats.n
    mean = stats.sx / n
    var = np.diag(stats.sxx) / n - mean**2
    sd = np.sqrt(np.maximum(var, 0.0))
    sd[sd == 0] = 1.0
    centered_gram = stats.sxx - n * np.outer(mean, mean)
    g_std = centered_gram / np.outer(sd, sd)
    b_std = (stats.sxy - mean * stats.sy) / sd
    return _StandardizedSystem(
        g_std=g_std, b_std=b_std, mean=mean, sd=sd, mean_y=stats.sy / n
    )


def _solve(system: _StandardizedSystem, alpha: float) -> np.ndarray:
    """Solve (G_std + alpha I) gamma_std = b_std for the standardized coefficients."""
    p = system.g_std.shape[0]
    return np.linalg.solve(system.g_std + alpha * np.eye(p), system.b_std)


def _heldout_r2(
    held: _ChromStats,
    gamma_std: np.ndarray,
    train_mean: np.ndarray,
    train_sd: np.ndarray,
    train_mean_y: float,
) -> float:
    """R^2 of the standardized model on a held-out chromosome, standardizing the
    held-out annotations with the TRAIN mean/sd.

    Prediction for raw annotations x is pred = train_mean_y + gamma_std . z, where
    z = (x - train_mean) / train_sd. Everything is expanded from the held-out
    chromosome's raw cross-product sufficient statistics; no per-SNP data needed.
    """
    n = held.n
    c = train_mean_y
    tm = train_mean
    ts = train_sd
    # SS over (y - c): sum (y - c)^2
    ss_res_y = held.syy - 2.0 * c * held.sy + n * c * c
    # sum z_i (y_i - c)  and  sum z_i z_i^T  in terms of raw stats
    z_r = (held.sxy - c * held.sx - tm * held.sy + n * c * tm) / ts
    zz = (
        held.sxx
        - np.outer(tm, held.sx)
        - np.outer(held.sx, tm)
        + n * np.outer(tm, tm)
    ) / np.outer(ts, ts)
    ss_res = (
        ss_res_y
        - 2.0 * float(gamma_std @ z_r)
        + float(gamma_std @ zz @ gamma_std)
    )
    mean_y_held = held.sy / n
    ss_tot = held.syy - n * mean_y_held * mean_y_held
    return 1.0 - ss_res / ss_tot


def _select_alpha_loco(
    per_chrom: dict[int, _ChromStats], alphas: tuple[float, ...]
) -> tuple[float, float, dict[str, float]]:
    """Pick alpha by leave-one-chromosome-out; return (alpha, mean_r2, r2_per_chrom)."""
    chroms = sorted(per_chrom)
    best_alpha = alphas[0]
    best_mean = -np.inf
    best_r2_per_chrom: dict[str, float] = {}
    for alpha in alphas:
        r2s: dict[str, float] = {}
        for held in chroms:
            train = _combine([per_chrom[c] for c in chroms if c != held])
            system = _standardized_system(train)
            gamma_std = _solve(system, alpha)
            r2s[str(held)] = _heldout_r2(
                per_chrom[held],
                gamma_std=gamma_std,
                train_mean=system.mean,
                train_sd=system.sd,
                train_mean_y=system.mean_y,
            )
        mean_r2 = float(np.mean(list(r2s.values())))
        if mean_r2 > best_mean:
            best_mean, best_alpha, best_r2_per_chrom = mean_r2, alpha, r2s
    return best_alpha, best_mean, best_r2_per_chrom


def _fit_all(
    per_chrom: dict[int, _ChromStats], alpha: float
) -> tuple[list[float], list[float], float]:
    """Fit on all chromosomes; return (gamma_raw, gamma_standardized, intercept)."""
    system = _standardized_system(_combine(list(per_chrom.values())))
    gamma_std = _solve(system, alpha)
    gamma_raw = gamma_std / system.sd
    intercept = float(system.mean_y - gamma_raw @ system.mean)
    return gamma_raw.tolist(), gamma_std.tolist(), intercept
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `pixi r pytest test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py -v`
Expected: PASS (`gamma_raw` within 1e-2 of truth, `mean_heldout_r2 > 0.999`, `n_variants == 800`).

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/build_system/task/annotation_weights/ridge_annotation_weights_task.py test_mecfs_bio/unit/build_system/task/annotation_weights/test_ridge_annotation_weights_task.py
git commit -m "feat: RidgeAnnotationWeightsTask (streaming Gram + LOCO alpha)"
```

---

### Task 4: Wire the real assets

**Files:**
- Create: `mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_annotations.py`
- Create: `mecfs_bio/assets/reference_data/polyfun/annotations/annotation_ridge_weights.py`
- Test: none new (wiring only; the Task logic is covered by Tasks 2–3). Import-smoke covered by `invoke green`'s import check.

**Interfaces:**
- Consumes: `DownloadFileTask`, `BuildBaselineLFAnnotationParquetTask`, `RidgeAnnotationWeightsTask`, and `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS` from `mecfs_bio/assets/reference_data/polyfun/precomputed_prior/polyfun_precomputed_prior.py`.
- Produces: `BASELINE_LF_ANNOTATION_TARBALL` (Task), `BASELINE_LF_ANNOTATION_MATRIX` (Task), `BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS` (Task).

- [ ] **Step 1: Wire the tarball + annotation-matrix tasks**

Create `baseline_lf_annotations.py`:

```python
"""Reference assets: the baseline-LF 2.2.UKB annotation tarball and the single
sorted annotation parquet built from it.

The tarball is ~11GB; it and the derived ~19M x 191 parquet are path_remap
candidates (large, few-file, rarely read) for relocation to /mnt/d.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.annotation_weights.build_baseline_lf_annotation_parquet_task import (
    BuildBaselineLFAnnotationParquetTask,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

BASELINE_LF_ANNOTATION_TARBALL = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="polyfun",
        sub_group="annotations",
        sub_folder=PurePath("raw"),
        id=AssetId("baseline_lf_2.2_ukb_annotations_tarball"),
        extension=".tar.gz",
    ),
    url="https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/baselineLF_v2.2.UKB.tar.gz",
    md5_hash=None,  # pin after the first real download (see Step 3)
)

BASELINE_LF_ANNOTATION_MATRIX = BuildBaselineLFAnnotationParquetTask.create(
    asset_id="baseline_lf_2.2_ukb_annotations",
    tarball_task=BASELINE_LF_ANNOTATION_TARBALL,
)
```

(If `ReferenceFileMeta` requires a `read_spec` for a non-dataframe file, check its signature and pass `read_spec=None` consistent with how other binary reference downloads are declared.)

- [ ] **Step 2: Wire the ridge-weights task**

Create `annotation_ridge_weights.py`:

```python
"""Reference asset: ridge surrogate weights (gamma_c) for the polyfun prior."""

from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.assets.reference_data.polyfun.precomputed_prior.polyfun_precomputed_prior import (
    COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    RidgeAnnotationWeightsTask,
)

BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS = RidgeAnnotationWeightsTask.create(
    asset_id="baseline_lf_2.2_ukb_annotation_ridge_weights",
    annotation_parquet_task=BASELINE_LF_ANNOTATION_MATRIX,
    snpvar_meta_task=COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
)
```

Note: `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS` produces columns `CHR, BP, SNP, A1, A2, snpvar_bin`; the ridge task selects only `SNP, snpvar_bin`, so it is compatible as-is.

- [ ] **Step 3: Run green (import + lint + type checks)**

Run: `pixi r invoke green 2>&1 | tee /tmp/green_annotation_weights.log`
Expected: PASS (imports resolve, types check, all Task 1–3 tests pass). Fix anything the log surfaces.

- [ ] **Step 4: Commit**

```bash
git add mecfs_bio/assets/reference_data/polyfun/annotations/baseline_lf_annotations.py mecfs_bio/assets/reference_data/polyfun/annotations/annotation_ridge_weights.py
git commit -m "feat: wire baseline-LF annotation + ridge-weights reference assets"
```

---

### Task 5: Real-data verification + md5 pin (manual, one-time)

**Not a CI test** — it downloads 11GB and processes ~19M variants. Run once locally to confirm the pipeline reproduces the spike, then pin the tarball md5.

- [ ] **Step 1: Build the annotation matrix and pin the md5**

Build the tarball + annotation matrix assets through the normal build entrypoint (the same mechanism other reference assets are materialized with — see how `COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS` is built). After the tarball lands, compute and pin its md5:

Run: `pixi r python -c "from mecfs_bio.util.download.verify import calc_md5_checksum; from pathlib import Path; print(calc_md5_checksum(Path('<path-to-downloaded-tarball>')))"`
Then set that value as `md5_hash=` in `baseline_lf_annotations.py` and commit (`chore: pin baseline-LF annotation tarball md5`).

- [ ] **Step 2: Build the ridge weights and check diagnostics**

Materialize `BASELINE_LF_ANNOTATION_RIDGE_WEIGHTS`, then inspect `diagnostics.json`:
- `mean_heldout_r2` should be ~0.85–0.88 (spike: 0.869).
- In `weights.parquet`, the largest `gamma_standardized` magnitudes should be conservation / coding / promoter-enhancer / MaxCPP annotations (spike top list). Record the check in `experiments/claude/polyfun_explain_probe/` for provenance.

- [ ] **Step 3: Note path_remap candidates**

Confirm with the user whether the 11GB tarball and the ~19M-row annotation parquet should be added to `default_runner_config.yaml` `path_remap` (relocated to /mnt/d). This config is gitignored and machine-local; do not commit it.

---

## Self-Review

**Spec coverage:**
- Component A (tarball download) → Task 4 Step 1 (+ md5 pin in Task 5). ✓
- Component B (annotation parquet, streaming sink, sorted, deduped, no sidecar; meta derived from tarball meta) → Task 2. ✓
- Component C (187→family map, 11 published-grounded families, keys==187 test, open_chromatin membership test) → Task 1. ✓
- Component D (streaming-Gram ridge, LOCO alpha, raw+standardized gamma, family col, diagnostics json, directory asset; meta derived from annotation-parquet meta) → Task 3. ✓
- Storage decision (single (CHR,BP)-sorted parquet) → Task 2. ✓
- Coefficient scale (gamma_raw + gamma_standardized, no mean_c/std_c) → Task 3 weights columns. ✓
- Meta-derivation pattern (create() reuses dep group/sub_group/sub_folder) → Task 2 & Task 3 create(). ✓
- Testing (injected tarball, injected fetch, synthetic linear recovery, family coverage) → Tasks 1–3. ✓
- Real-data reproduce-the-spike gate → Task 5. ✓

**Placeholder scan:** No "TBD/TODO/handle edge cases". The one deferred value (`md5_hash=None`) is an explicit, commanded pin step (Task 5 Step 1), not a vague placeholder. The ridge scoring helpers are a single clean held-out-R^2 path (`_heldout_r2`) with no dead code.

**Type consistency:** `BuildBaselineLFAnnotationParquetTask.create(asset_id, tarball_task)`, `RidgeAnnotationWeightsTask.create(asset_id, annotation_parquet_task, snpvar_meta_task, alphas)`, and the constants `ANNOT_KEY_COLUMNS`, `WEIGHTS_PARQUET_FILENAME`, `DIAGNOSTICS_JSON_FILENAME`, `ANNOTATION_COL`, `GAMMA_RAW_COL`, `SNPVAR_COL` are referenced consistently across the task modules, the wiring, and the tests. Both `create()`s derive their meta from the dependency's `ReferenceFileMeta` (reusing group/sub_group/sub_folder) and `raise` on an unknown meta. `_select_alpha_loco(per_chrom, alphas)` and `_heldout_r2(...)` signatures match their call sites; every helper (incl. `_StandardizedSystem`) and test helper carries parameter+return annotations. Annotation columns are always "all parquet columns except `ANNOT_KEY_COLUMNS`".

**Family taxonomy provenance:** 11 families — 7 names verbatim from Weissbrod 2020's Supplementary Note grouping (non_synonymous, coding, conserved, promoter_or_enhancer, histone_marks, repressed, other), maf_bins + ld_related_continuous from Gazal 2017, molecular_qtl from Hormozdiari 2018, and open_chromatin as the single documented deviation (DHS/FetalDHS/DGF split out of polyfun's "others"). Cited in the module docstring.

**Known risk to watch during implementation:** the LOCO held-out R^2 sufficient-statistics algebra in `_r2_from_stats` is the subtlest piece — Task 3's synthetic test (noise-free exact-linear → `mean_heldout_r2 > 0.999`) is the guard. If that assertion is off, verify `_r2_from_stats` against a brute-force per-row computation on the small test data before touching anything else.
