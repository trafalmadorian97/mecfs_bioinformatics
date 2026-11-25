# gnomAD and Pan-UKBB Allele-Frequency Panels Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let GenomeReferenceHarmonizationTask resolve palindromes and ambiguous indels against
a chosen ancestry of a multi-ancestry allele-frequency panel, and add three such panels: gnomAD
v2.1.1 genomes (hg19), gnomAD v4.1 genomes (hg38, defined but built on request) and the Pan-UKBB
variant manifest (hg19).

**Architecture:** Panels declare their (ancestry, AF column) pairs in
HarmonizableReferenceTableMeta, and each harmonization names its ancestry; the loader selects
that column as AF. Both new panel sources stream large sorted tables through one shared batch
checker (contig mapping, ACGT alleles, order, uniqueness, REF vs FASTA) into parquet. gnomAD is
built one chromosome per Task (bcftools query over HTTPS, resumable), then concatenated;
Pan-UKBB is one Task over a pinned download.

**Tech Stack:** Python 3, polars (scan_csv + collect_batches for reading), pyarrow 25 (ParquetWriter), numpy
memmap FASTA gather, bcftools 1.24 / htslib (pixi), attrs, pytest, the repo's build system.

**Spec:** experiments/claude/design_specs/2026-10-02-gnomad-allele-frequency-panel-design.md

## Global Constraints

- Run everything through pixi: `pixi r python ...`, `pixi r python -m pytest ...`, `pixi r invoke green`.
- After each code task: `pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_<task>.log`, then read the log; the exit code alone proves nothing under testmon.
- Docstrings: no backticks around inline code, no RST.
- Imports at the top of the file, never inside functions or tests.
- No mocks or monkeypatching; inject dependencies. Tests are Task-level. Never match on error text (`pytest.raises(AssertionError)` with no `match=`).
- Task classes and reusable machinery live under mecfs_bio/build_system; URLs, release constants and asset instances live under mecfs_bio/assets; shared Literals live in mecfs_bio/constants. Never hardcode a URL in build_system.
- Prefer free helper functions over methods; frozen attrs classes (not tuples) for multi-value returns; Literal types for enumerable strings; assert preconditions with clear messages.
- Use the column constants GWASLAB_CHROM_COL ("CHR") and GWASLAB_POS_COL ("POS") from mecfs_bio/constants/gwaslab_constants.py and PANEL_REF_COL / PANEL_ALT_COL / PANEL_AF_COL from reference_panel_task.py.
- Existing harmonized outputs must not change (1000 Genomes panels pass panel_ancestry="eur" and keep their "AF" column).
- Do not wrap the new panel Tasks in DiscardDepsWrapper.
- Parquet encoding for panels: AF columns Float32 with dictionary encoding; gnomAD AN columns Int32 with byte-stream-split; zstd at the default level.
- Pan-UKBB: expected REF/FASTA mismatches = 295 exactly; zero FASTA-ambiguous rows; manifest md5 e70ebc8289f762dd8d5086f54e766654.
- gnomAD: zero REF/FASTA mismatches over pure A/C/G/T spans; records over N/IUPAC are dropped and counted.
- Comments describe current code, not history.

## Review Focus

1. Split multiallelic sites: two gnomAD records at one position with different ALT (A>G, A>T) must both be kept, not flagged as duplicates, including when the pair straddles a batch boundary. Pinned in Task 4's happy-path test, run with default and one-row batches.
2. Long REF alleles (gnomAD has deletions of hundreds of bases): the FASTA gather must classify them, not crash or allocate unboundedly. Pinned in Task 4's happy-path test with a 30-base deletion.
3. A REF span running past the end of its FASTA contig (wrong build or contig naming): must fail loudly, not be silently counted as FASTA-ambiguous. Pinned in Task 3 (test_ref_span_beyond_contig_end_fails).
4. A chromosome with no PASS record polymorphic in a main group (bcftools writes an empty TSV): must fail with a clear message rather than a CSV-reader error or an empty panel. Pinned in Task 4 (test_chromosome_without_polymorphic_pass_records_fails).
5. A panel whose chosen ancestry column is null at a site (AN = 0 in that group): the variant must be treated as absent from the panel, not resolved with AF 0. Pinned in Task 1 (test_null_chosen_ancestry_frequency_counts_as_absent).

---

## File Structure

Create:
- `mecfs_bio/constants/allele_frequency_panel_constants.py`: GnomadGroup, PanUkbbGroup, PanelAncestry Literals; panel_af_col / panel_an_col.
- `mecfs_bio/build_system/meta/reference_meta/panel_allele_frequency_columns.py`: PanelAncestryColumn, PanelAlleleFrequencyColumns.
- `mecfs_bio/build_system/task/genome_reference_harmonization/panel_batch_checks.py`: shared batch checks, write_checked_panel, fetch_file_path.
- `mecfs_bio/build_system/task/genome_reference_harmonization/pan_ukbb/__init__.py` (empty)
- `mecfs_bio/build_system/task/genome_reference_harmonization/pan_ukbb/pan_ukbb_allele_frequency_panel_task.py`
- `mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/__init__.py` (empty)
- `mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/gnomad_release.py`
- `mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/gnomad_chromosome_allele_frequency_task.py`
- `mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/gnomad_allele_frequency_panel_task.py`
- `mecfs_bio/build_system/task_generator/gnomad_allele_frequency_panel_task_generator.py` (builds the parts and the panel together)
- `mecfs_bio/assets/reference_data/gnomad/__init__.py` (empty), `gnomad_releases.py`, `gnomad_allele_frequency_panels.py`
- `mecfs_bio/assets/reference_data/pan_ukbb/__init__.py` (empty), `pan_ukbb_variant_manifest.py`, `pan_ukbb_allele_frequencies.py`
- Tests in `test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/`: `test_pan_ukbb_allele_frequency_panel_task.py`, `test_gnomad_chromosome_allele_frequency_task.py`, `test_gnomad_allele_frequency_panel_task.py`
- Experiments in `experiments/claude/gnomad_af_reference/`: `check_release_headers.py`, `byte_identity_1000g.py`, `build_panel.py`, `decode_me_panel_comparison.py`

Modify:
- `mecfs_bio/build_system/meta/reference_meta/harmonizable_reference_table_meta.py`: add allele_frequency_columns.
- `mecfs_bio/build_system/task/genome_reference_harmonization/reference_panel_task.py`: required ancestry in create; pipe comment.
- `mecfs_bio/build_system/task/genome_reference_harmonization/genome_reference_harmonization_task.py`: panel_ancestry, panel_af_col, PanelTable, ParquetPanelLoader.af_col.
- `mecfs_bio/build_system/task/genome_reference_harmonization/fasta.py`: gwaslab_code_to_contig_name, reference_is_acgt, shared gather.
- `mecfs_bio/build_system/task/dataframe_output.py`: ParquetEncoding, parquet_encoding, open_parquet_writer.
- Call sites: `mecfs_bio/assets/reference_data/thousand_genomes/eur_panel_allele_frequencies.py`, `mecfs_bio/asset_generator/annovar_37_basic_rsid_assignment.py`, `mecfs_bio/assets/gwas/inflammatory_bowel_disease/liu_et_al_2023/processed_gwas_data/liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized.py`, `test_mecfs_bio/system/test_harmonize_drop_ambiguous.py`.
- Experiments using the changed helpers: `experiments/claude/genome_reference_harmonization/{survey_trust,v3_other_truth_sets,v3_suspicious_calibration,tune_ambiguous_indel_rules}.py`.
- Tests: `genome_reference_fixtures.py`, `test_genome_reference_harmonization_task.py`, `test_reference_panel_task.py`, `test_mecfs_bio/unit/build_system/task/test_dataframe_output.py`.

---

### Task 1: Panel ancestry declaration and harmonizer column selection

**Files:**
- Create: `mecfs_bio/constants/allele_frequency_panel_constants.py`
- Create: `mecfs_bio/build_system/meta/reference_meta/panel_allele_frequency_columns.py`
- Modify: `mecfs_bio/build_system/meta/reference_meta/harmonizable_reference_table_meta.py`
- Modify: `mecfs_bio/build_system/task/genome_reference_harmonization/reference_panel_task.py`
- Modify: `mecfs_bio/build_system/task/genome_reference_harmonization/genome_reference_harmonization_task.py`
- Modify: the four call sites and four experiments listed above
- Test: `test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/genome_reference_fixtures.py`, `test_genome_reference_harmonization_task.py`, `test_reference_panel_task.py`

**Interfaces:**
- Produces:
  - `GnomadGroup`, `PanUkbbGroup`, `PanelAncestry` (Literals), `panel_af_col(ancestry: PanelAncestry) -> str` ("AF_" + ancestry), `panel_an_col(ancestry: PanelAncestry) -> str` ("AN_" + ancestry) in `mecfs_bio.constants.allele_frequency_panel_constants`.
  - `PanelAncestryColumn(ancestry, column)`, `PanelAlleleFrequencyColumns(entries)` with `.ancestries -> list[PanelAncestry]`, `.column_for(ancestry) -> str`, classmethods `.single(ancestry, column)` and `.prefixed(ancestries: Sequence[PanelAncestry])`.
  - `HarmonizableReferenceTableMeta.allele_frequency_columns: PanelAlleleFrequencyColumns | None = None`.
  - `ReferencePanelAlleleFrequencyTask.create(vcf_task, asset_id, build, ancestry: PanelAncestry)`.
  - `GenomeReferenceHarmonizationTask.create(asset_id, sumstats_task, fasta_task, panel_task, panel_ancestry: PanelAncestry, options=..., pipe=...)`; Task field `panel_af_col: str`.
  - `PanelTable(path: Path, af_col: str)`; `ParquetPanelLoader(panel_path, af_col, chrom)`; `count_trust_evidence_genome_wide(sumstats, chromosomes, fasta, panel: PanelTable, options)`; `resolve_chromosome_rows(sumstats, context, panel: PanelTable)`.
  - Fixtures: `write_fasta(directory, sequences=SYNTHETIC_SEQUENCES)`, `AncestryPanelRecord`, `EUR_PANEL_COLUMNS`, `ANCESTRY_PANEL_COLUMNS`, `harmonization_task(panel_columns, panel_ancestry, options=TEST_OPTIONS, pipe=IdentityPipe())`, `run_harmonization(..., ancestry_panel=None, panel_ancestry="eur")`.

- [ ] **Step 1: Create the vocabulary module**

`mecfs_bio/constants/allele_frequency_panel_constants.py`:

```python
"""
Ancestry vocabularies of the allele-frequency panels used by genome-reference harmonization.

A panel declares which ancestries it carries, and each harmonization names the one that
matches its GWAS. The vocabularies are deliberately disjoint: "eur" is the 1000 Genomes EUR
super-population, gnomAD groups are bare gnomAD codes (nfe, afr, ...), and Pan-UKBB groups
carry a ukb_ prefix. "eur", "nfe" and "ukb_eur" are different populations and never synonyms.

Panels built from multi-ancestry sources name their columns AF_<ancestry> and
AN_<ancestry>.
"""

from typing import Literal

GnomadGroup = Literal[
    "afr",
    "ami",
    "amr",
    "asj",
    "eas",
    "fin",
    "mid",
    "nfe",
    "nfe_est",
    "nfe_nwe",
    "nfe_onf",
    "nfe_seu",
    "oth",
    "remaining",
    "sas",
]
PanUkbbGroup = Literal[
    "ukb_afr", "ukb_amr", "ukb_csa", "ukb_eas", "ukb_eur", "ukb_mid"
]
ThousandGenomesSuperPopulation = Literal["eur"]
PanelAncestry = GnomadGroup | PanUkbbGroup | ThousandGenomesSuperPopulation

PANEL_AF_PREFIX = "AF_"
PANEL_AN_PREFIX = "AN_"


def panel_af_col(ancestry: PanelAncestry) -> str:
    """Name of a multi-ancestry panel's allele-frequency column for one ancestry."""
    return PANEL_AF_PREFIX + ancestry


def panel_an_col(ancestry: PanelAncestry) -> str:
    """Name of a multi-ancestry panel's allele-number column for one ancestry."""
    return PANEL_AN_PREFIX + ancestry
```

- [ ] **Step 2: Create PanelAlleleFrequencyColumns**

`mecfs_bio/build_system/meta/reference_meta/panel_allele_frequency_columns.py`:

```python
"""Which allele-frequency column of a reference panel table holds which ancestry."""

from collections.abc import Sequence

from attrs import frozen

from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanelAncestry,
    panel_af_col,
)


@frozen(slots=True)
class PanelAncestryColumn:
    ancestry: PanelAncestry
    column: str


@frozen(slots=True)
class PanelAlleleFrequencyColumns:
    """The ancestries a panel carries, each with its allele-frequency column.

    Entries are a tuple so the declaration stays hashable on a frozen Task's meta.
    """

    entries: tuple[PanelAncestryColumn, ...]

    def __attrs_post_init__(self):
        ancestries = self.ancestries
        assert ancestries, "a panel must declare at least one ancestry"
        assert len(set(ancestries)) == len(ancestries), (
            f"duplicate panel ancestries: {ancestries}"
        )

    @property
    def ancestries(self) -> list[PanelAncestry]:
        return [entry.ancestry for entry in self.entries]

    def column_for(self, ancestry: PanelAncestry) -> str:
        matches = [entry.column for entry in self.entries if entry.ancestry == ancestry]
        assert matches, (
            f"the panel has no allele frequency for ancestry {ancestry!r}; "
            f"available: {self.ancestries}"
        )
        return matches[0]

    @classmethod
    def single(
        cls, ancestry: PanelAncestry, column: str
    ) -> "PanelAlleleFrequencyColumns":
        return cls((PanelAncestryColumn(ancestry=ancestry, column=column),))

    @classmethod
    def prefixed(
        cls, ancestries: Sequence[PanelAncestry]
    ) -> "PanelAlleleFrequencyColumns":
        """One AF_<ancestry> column per ancestry, as multi-ancestry panels name them."""
        return cls(
            tuple(
                PanelAncestryColumn(ancestry=ancestry, column=panel_af_col(ancestry))
                for ancestry in ancestries
            )
        )
```

- [ ] **Step 3: Add the meta field**

In `harmonizable_reference_table_meta.py`, import `PanelAlleleFrequencyColumns` from `mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns` and add after `harmonization_info`:

```python
    harmonization_info: HarmonizationInfo | None = None
    # Set by allele-frequency panels: which column holds which ancestry's frequency.
    allele_frequency_columns: PanelAlleleFrequencyColumns | None = None
```

- [ ] **Step 4: Update the fixtures so the harmonization tests can choose a panel ancestry**

In `genome_reference_fixtures.py`:

1. Add imports (top of file):

```python
from collections.abc import Mapping, Sequence

from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)
from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanelAncestry,
    panel_af_col,
)
```

2. Rename `_write_fasta` to a public, parametrized `write_fasta`, keeping the old sequences as the default:

```python
SYNTHETIC_SEQUENCES: Mapping[str, str] = {"chr1": CHR1_SEQUENCE, "chr2": CHR2_SEQUENCE}


def write_fasta(
    directory: Path, sequences: Mapping[str, str] = SYNTHETIC_SEQUENCES
) -> None:
    """Write an uncompressed, faidx-indexed FASTA with LINE_WIDTH-base lines."""
    directory.mkdir()
    fasta_path = directory / FASTA_FILENAME
    lines: list[str] = []
    for name, sequence in sequences.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[start : start + LINE_WIDTH]
            for start in range(0, len(sequence), LINE_WIDTH)
        )
    fasta_path.write_text("\n".join(lines) + "\n")
    pysam.faidx(str(fasta_path))
```

3. Add the two-ancestry panel record, column declarations and writer:

```python
EUR_PANEL_COLUMNS = PanelAlleleFrequencyColumns.single("eur", PANEL_AF_COL)
ANCESTRY_PANEL_COLUMNS = PanelAlleleFrequencyColumns.prefixed(("nfe", "nfe_nwe"))


@frozen(slots=True)
class AncestryPanelRecord:
    """A panel record with gnomAD-style AF_nfe and AF_nfe_nwe columns (None = AN 0)."""

    pos: int
    ref: str
    alt: str
    nfe: float | None
    nfe_nwe: float | None
    chrom: int = 1


def _write_ancestry_panel(path: Path, records: Sequence[AncestryPanelRecord]) -> None:
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [r.chrom for r in records],
            GWASLAB_POS_COL: [r.pos for r in records],
            PANEL_REF_COL: [r.ref for r in records],
            PANEL_ALT_COL: [r.alt for r in records],
            panel_af_col("nfe"): [r.nfe for r in records],
            panel_af_col("nfe_nwe"): [r.nfe_nwe for r in records],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            panel_af_col("nfe"): pl.Float32,
            panel_af_col("nfe_nwe"): pl.Float32,
        },
    ).sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL).write_parquet(path)
```

4. Extract the Task construction into a public builder and use it from run_harmonization:

```python
def harmonization_task(
    panel_columns: PanelAlleleFrequencyColumns | None,
    panel_ancestry: PanelAncestry,
    options: GenomeReferenceHarmonizationOptions = TEST_OPTIONS,
    pipe: DataProcessingPipe = IdentityPipe(),
) -> GenomeReferenceHarmonizationTask:
    """The Task under test, wired to FakeTask inputs whose ids run_harmonization serves."""
    parquet_spec = DataFrameReadSpec(DataFrameParquetFormat())
    return GenomeReferenceHarmonizationTask.create(
        asset_id="harmonized",
        sumstats_task=FakeTask(
            FilteredGWASDataMeta(
                id=AssetId(_SUMSTATS_ID),
                trait="synthetic_trait",
                project="synthetic_project",
                sub_dir="processed",
                read_spec=parquet_spec,
            )
        ),
        fasta_task=FakeTask(
            FASTAMeta(
                group="genome_sequence",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                id=AssetId(_FASTA_ID),
                build="19",
            )
        ),
        panel_task=FakeTask(
            HarmonizableReferenceTableMeta(
                group="reference_panel_allele_frequencies",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                extension=".parquet",
                id=AssetId(_PANEL_ID),
                read_spec=parquet_spec,
                harmonization_info=HarmonizationInfo(
                    build="19", ref_allele_col=PANEL_REF_COL, pos_col=GWASLAB_POS_COL
                ),
                allele_frequency_columns=panel_columns,
            )
        ),
        panel_ancestry=panel_ancestry,
        options=options,
        pipe=pipe,
    )


def run_harmonization(
    work_dir: Path,
    sumstats: pl.DataFrame,
    panel: Sequence[PanelRecord] = (),
    options: GenomeReferenceHarmonizationOptions = TEST_OPTIONS,
    pipe: DataProcessingPipe = IdentityPipe(),
    ancestry_panel: Sequence[AncestryPanelRecord] | None = None,
    panel_ancestry: PanelAncestry = "eur",
) -> pl.DataFrame:
    """Execute the Task on synthetic inputs in a fresh work_dir and return the output table.

    By default the panel is the single-ancestry "eur" panel built from panel. Pass
    ancestry_panel to use a panel with AF_nfe and AF_nfe_nwe columns instead.
    """
    work_dir.mkdir(parents=True)
    sumstats_path = work_dir / "sumstats.parquet"
    sumstats.write_parquet(sumstats_path)
    fasta_dir = work_dir / "fasta"
    write_fasta(fasta_dir)
    panel_path = work_dir / "panel.parquet"
    if ancestry_panel is None:
        _write_panel(panel_path, panel)
        panel_columns = EUR_PANEL_COLUMNS
    else:
        _write_ancestry_panel(panel_path, ancestry_panel)
        panel_columns = ANCESTRY_PANEL_COLUMNS
    task = harmonization_task(
        panel_columns=panel_columns,
        panel_ancestry=panel_ancestry,
        options=options,
        pipe=pipe,
    )
    assets: dict[str, Asset] = {
        _SUMSTATS_ID: FileAsset(sumstats_path),
        _FASTA_ID: DirectoryAsset(fasta_dir),
        _PANEL_ID: FileAsset(panel_path),
    }

    def fetch(asset_id: AssetId) -> Asset:
        return assets[asset_id]

    scratch = work_dir / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)
```

Remove the old `parquet_spec` / inline `GenomeReferenceHarmonizationTask.create` block from run_harmonization (now in harmonization_task).

- [ ] **Step 5: Write the failing tests**

Append to `test_genome_reference_harmonization_task.py` (after `_STRINGENT_OPTIONS` is defined, i.e. at the end of the file). Add `AncestryPanelRecord`, `EUR_PANEL_COLUMNS`, `harmonization_task` to the fixture import list and `from mecfs_bio.constants.allele_frequency_panel_constants import PanelAncestry` to the imports.

```python
_ANCESTRY_PALINDROME = Variant(pos=13, ea="A", nea="T", eaf=0.1, beta=0.3)
# nfe is on the same side of 0.5 as the EAF (keep); nfe_nwe is on the other (strand flip).
_ANCESTRY_PALINDROME_PANEL = [
    AncestryPanelRecord(pos=13, ref="T", alt="A", nfe=0.15, nfe_nwe=0.9)
]
_ANCESTRY_CASES: list[tuple[PanelAncestry, float, float]] = [
    ("nfe", 0.3, 0.1),
    ("nfe_nwe", -0.3, 0.9),
]


@pytest.mark.parametrize("ancestry, expected_beta, expected_eaf", _ANCESTRY_CASES)
def test_palindrome_is_resolved_with_the_chosen_ancestry_column(
    tmp_path: Path, ancestry: PanelAncestry, expected_beta: float, expected_eaf: float
) -> None:
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV, _ANCESTRY_PALINDROME]),
        ancestry_panel=_ANCESTRY_PALINDROME_PANEL,
        panel_ancestry=ancestry,
    )
    row = row_at(result, 13)
    assert (row[EA], row[NEA]) == ("A", "T")
    assert (row[BETA], row[EAF]) == (
        pytest.approx(expected_beta),
        pytest.approx(expected_eaf),
    )


def test_null_chosen_ancestry_frequency_counts_as_absent(tmp_path: Path) -> None:
    # The nfe reading fits (see test_untrusted_ambiguous_indels_follow_the_stringent_rules,
    # pos 33); nfe_nwe has no frequency (AN 0), so under it the indel is not in the panel.
    variant = Variant(pos=33, ea="T", nea="TT", eaf=0.3)
    panel = [AncestryPanelRecord(pos=33, ref="TT", alt="T", nfe=0.32, nfe_nwe=None)]
    ancestries: list[PanelAncestry] = ["nfe", "nfe_nwe"]
    kept = {
        ancestry: 33
        in positions(
            run_harmonization(
                tmp_path / ancestry,
                sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV, variant]),
                ancestry_panel=panel,
                panel_ancestry=ancestry,
                options=_STRINGENT_OPTIONS,
            )
        )
        for ancestry in ancestries
    }
    assert kept == {"nfe": True, "nfe_nwe": False}


def test_undeclared_panel_ancestry_is_rejected() -> None:
    with pytest.raises(AssertionError):
        harmonization_task(panel_columns=EUR_PANEL_COLUMNS, panel_ancestry="nfe")


def test_panel_without_declared_columns_is_rejected() -> None:
    with pytest.raises(AssertionError):
        harmonization_task(panel_columns=None, panel_ancestry="eur")
```

- [ ] **Step 6: Run the tests to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py -v`
Expected: every test fails (TypeError: create() got an unexpected keyword argument 'panel_ancestry'), since the fixture now passes it.

- [ ] **Step 7: Implement the harmonizer changes**

In `genome_reference_harmonization_task.py`:

1. Imports: add `from mecfs_bio.constants.allele_frequency_panel_constants import PanelAncestry`.
2. Module docstring: append the sentence "The panel's allele-frequency column is chosen per harmonization through panel_ancestry, from the ancestries the panel's meta declares."
3. Replace ParquetPanelLoader and add PanelTable:

```python
@frozen(slots=True)
class PanelTable:
    """A reference panel parquet and the allele-frequency column this harmonization reads."""

    path: Path
    af_col: str


@frozen(slots=True)
class ParquetPanelLoader:
    """Reads one chromosome's panel rows at requested positions from the panel parquet.

    The chosen af_col is returned as PANEL_AF_COL. Rows where it is null (nobody in that
    ancestry was called at the site) are dropped, so the record counts as absent.
    """

    panel_path: Path
    af_col: str
    chrom: int

    def __call__(self, positions: pl.Series) -> pl.DataFrame:
        wanted = positions.cast(pl.Int32).unique().to_frame(GWASLAB_POS_COL).lazy()
        return (
            pl.scan_parquet(self.panel_path)
            .filter(pl.col(GWASLAB_CHROM_COL) == self.chrom)
            .join(wanted, on=GWASLAB_POS_COL, how="semi")
            .select(
                pl.col(GWASLAB_POS_COL).cast(pl.Int64),
                PANEL_REF_COL,
                PANEL_ALT_COL,
                pl.col(self.af_col).alias(PANEL_AF_COL),
            )
            .filter(pl.col(PANEL_AF_COL).is_not_null())
            .collect()
        )
```

4. `count_trust_evidence_genome_wide`: replace parameter `panel_path: Path` by `panel: PanelTable`, and the loader construction by `ParquetPanelLoader(panel_path=panel.path, af_col=panel.af_col, chrom=chrom)`.
5. `resolve_chromosome_rows(sumstats, context, panel: PanelTable)`: `load_panel=ParquetPanelLoader(panel_path=panel.path, af_col=panel.af_col, chrom=context.chrom)`.
6. `_write_chromosome_part(sumstats, context, panel: PanelTable, parts_dir)`: pass `panel` through.
7. Add the resolver next to resolve_harmonized_build:

```python
def resolve_panel_af_col(panel_task: Task, panel_ancestry: PanelAncestry) -> str:
    """The panel column holding panel_ancestry's allele frequency, from the panel's meta."""
    panel_meta = panel_task.meta
    assert (
        isinstance(panel_meta, HarmonizableReferenceTableMeta)
        and panel_meta.allele_frequency_columns is not None
    ), "panel_task must carry HarmonizableReferenceTableMeta with allele_frequency_columns"
    return panel_meta.allele_frequency_columns.column_for(panel_ancestry)
```

8. Task class: add the field `panel_af_col: str` directly after `panel_task: Task`. In execute, replace `panel_path = load_panel_path(fetch, self.panel_task)` with

```python
        panel = PanelTable(
            path=load_panel_path(fetch, self.panel_task), af_col=self.panel_af_col
        )
```

and pass `panel` (instead of `panel_path`) to `count_trust_evidence_genome_wide` and `_write_chromosome_part(..., panel=panel, parts_dir=parts_dir)`.

9. `create`: add the required parameter `panel_ancestry: PanelAncestry` after `panel_task: Task` (before `options`), and pass `panel_af_col=resolve_panel_af_col(panel_task, panel_ancestry)` to `cls(...)`. Add to the create docstring (write one if absent): "panel_ancestry is required, never defaulted: a silent default is how the wrong ancestry slips in."

- [ ] **Step 8: Declare the 1000 Genomes panel ancestry and comment on the pipe**

In `reference_panel_task.py`:

1. Imports: `PanelAlleleFrequencyColumns` and `PanelAncestry`.
2. `create(cls, vcf_task: Task, asset_id: str, build: GenomeBuild, ancestry: PanelAncestry)`; in the meta add `allele_frequency_columns=PanelAlleleFrequencyColumns.single(ancestry, PANEL_AF_COL)`.
3. Above the `execute_command` call in `_write_sites_tsv`, add:

```python
    # execute_command runs this through sh, where a pipeline's exit status is the last
    # command's: if bcftools view failed partway (for example on a truncated VCF), query
    # could still exit 0 after writing a truncated TSV. A single "bcftools query -i"
    # process, as GnomadChromosomeAlleleFrequencyTask uses, would fail instead.
```

In `test_reference_panel_task.py`, pass `ancestry="eur"` to `ReferencePanelAlleleFrequencyTask.create`.

- [ ] **Step 9: Update the call sites**

- `mecfs_bio/assets/reference_data/thousand_genomes/eur_panel_allele_frequencies.py`: add `ancestry="eur"` to both `ReferencePanelAlleleFrequencyTask.create` calls.
- `mecfs_bio/asset_generator/annovar_37_basic_rsid_assignment.py`, `.../liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized.py`, `test_mecfs_bio/system/test_harmonize_drop_ambiguous.py`: add `panel_ancestry="eur",` after `panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,`.
- Experiments (not type-checked, so update them by hand): in `experiments/claude/genome_reference_harmonization/survey_trust.py`, `v3_other_truth_sets.py`, `v3_suspicious_calibration.py`, `tune_ambiguous_indel_rules.py`, wherever a panel path is passed to `count_trust_evidence_genome_wide` or `resolve_chromosome_rows`, pass `PanelTable(path=<that path>, af_col=PANEL_AF_COL)` instead, and add `af_col=PANEL_AF_COL` to any `ParquetPanelLoader(...)`; also add `panel_ancestry="eur"` to the `GenomeReferenceHarmonizationTask.create` call in `measure_harmonization_memory.py`. Import `PanelTable` from the harmonization task module and `PANEL_AF_COL` from reference_panel_task. Leave `compare_with_gwaslab.py` alone: its docstring records it no longer imports. Check each edited experiment still imports: `pixi r python -c "import experiments.claude.genome_reference_harmonization.survey_trust"` (repeat per module).

Run `grep -rn "GenomeReferenceHarmonizationTask.create\|ReferencePanelAlleleFrequencyTask.create\|ParquetPanelLoader(\|count_trust_evidence_genome_wide(\|resolve_chromosome_rows(" --include=*.py mecfs_bio test_mecfs_bio experiments` and confirm every hit is updated (except compare_with_gwaslab.py).

- [ ] **Step 10: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/ -v`
Expected: all PASS, including the four new tests and every pre-existing test unchanged.

- [ ] **Step 11: Run green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_task1.log
git add -A mecfs_bio test_mecfs_bio experiments/claude/genome_reference_harmonization
git commit -m "Let harmonization choose a panel ancestry column

Panels declare (ancestry, AF column) pairs in HarmonizableReferenceTableMeta;
GenomeReferenceHarmonizationTask.create takes a required panel_ancestry and
reads that column, treating null frequencies as absent. 1000 Genomes panels
declare eur -> AF, so existing outputs are unchanged.

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: One parquet encoding helper for table and streaming writers

**Files:**
- Modify: `mecfs_bio/build_system/task/dataframe_output.py`
- Test: `test_mecfs_bio/unit/build_system/task/test_dataframe_output.py`

**Interfaces:**
- Produces: `ParquetEncoding(compression, compression_level, use_byte_stream_split: list[str] | bool, use_dictionary: list[str] | bool)`; `parquet_encoding(column_names: Sequence[str], compression: ParquetCompression, compression_level: int | None, byte_stream_split_columns: Sequence[str]) -> ParquetEncoding`; `open_parquet_writer(out_path: Path, schema: pyarrow.Schema, encoding: ParquetEncoding) -> pyarrow.parquet.ParquetWriter`. write_parquet_table keeps its signature.

- [ ] **Step 1: Write the failing test**

Append to `test_dataframe_output.py` (add `import pyarrow as pa` and `open_parquet_writer`, `parquet_encoding` to the imports):

```python
def test_streaming_writer_splits_named_int_columns_and_keeps_dictionary_elsewhere(
    tmp_path: Path,
):
    table = pa.table(
        {
            "count": pa.array([1, 2, 3, 4], pa.int32()),
            "label": ["a", "b", "a", "b"],
        }
    )
    encoding = parquet_encoding(
        table.schema.names,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=["count"],
    )
    path = tmp_path / "out.parquet"
    with open_parquet_writer(path, table.schema, encoding) as writer:
        writer.write_table(table)
        writer.write_table(table)
    assert "BYTE_STREAM_SPLIT" in _encodings(path, "count")
    assert "RLE_DICTIONARY" in _encodings(path, "label")
    assert pq.read_table(path).num_rows == 8
```

- [ ] **Step 2: Run it to verify it fails**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/test_dataframe_output.py -v`
Expected: FAIL with ImportError (cannot import name 'open_parquet_writer').

- [ ] **Step 3: Implement**

In `dataframe_output.py`, add after `OutFormat`:

```python
@frozen(slots=True)
class ParquetEncoding:
    """pyarrow parquet writer encoding arguments, as derived by parquet_encoding."""

    compression: ParquetCompression
    compression_level: int | None
    use_byte_stream_split: list[str] | bool
    use_dictionary: list[str] | bool


def parquet_encoding(
    column_names: Sequence[str],
    compression: ParquetCompression,
    compression_level: int | None,
    byte_stream_split_columns: Sequence[str],
) -> ParquetEncoding:
    """Writer encoding with BYTE_STREAM_SPLIT on the named columns and dictionary elsewhere.

    Dictionary encoding takes precedence over BYTE_STREAM_SPLIT in the parquet writer: a
    column left dictionary-enabled is written as RLE_DICTIONARY and the requested split is
    silently dropped. Dictionary encoding is therefore disabled on exactly the split
    columns and left on for the rest, where it is what makes low-cardinality columns small.
    The split applies to any fixed-width column (pyarrow 25 supports integer as well as
    floating-point types).

    Pass an empty byte_stream_split_columns to disable the split entirely.
    """
    split_columns = list(byte_stream_split_columns)
    missing = set(split_columns) - set(column_names)
    assert not missing, f"byte_stream_split_columns not in frame: {missing}"
    other_columns = [name for name in column_names if name not in set(split_columns)]
    return ParquetEncoding(
        compression=compression,
        compression_level=compression_level,
        use_byte_stream_split=split_columns if split_columns else False,
        use_dictionary=other_columns if split_columns else True,
    )


def open_parquet_writer(
    out_path: Path, schema: pyarrow.Schema, encoding: ParquetEncoding
) -> pyarrow.parquet.ParquetWriter:
    """A streaming parquet writer with the given encoding; close it (or use with)."""
    return pyarrow.parquet.ParquetWriter(
        str(out_path),
        schema,
        compression=encoding.compression,
        compression_level=encoding.compression_level,
        use_byte_stream_split=encoding.use_byte_stream_split,
        use_dictionary=encoding.use_dictionary,
    )
```

Replace the body of `write_parquet_table` (keep its signature) with:

```python
    """Write an arrow table to parquet with explicit encoding control (see parquet_encoding)."""
    encoding = parquet_encoding(
        table.schema.names,
        compression=compression,
        compression_level=compression_level,
        byte_stream_split_columns=byte_stream_split_columns,
    )
    pyarrow.parquet.write_table(
        table,
        out_path,
        compression=encoding.compression,
        compression_level=encoding.compression_level,
        use_byte_stream_split=encoding.use_byte_stream_split,
        use_dictionary=encoding.use_dictionary,
    )
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/test_dataframe_output.py -v`
Expected: all PASS (new test and the existing split tests).

- [ ] **Step 5: Run green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_task2.log
git add mecfs_bio/build_system/task/dataframe_output.py test_mecfs_bio/unit/build_system/task/test_dataframe_output.py
git commit -m "Share parquet encoding between table and streaming writers

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Shared panel batch checks and the Pan-UKBB panel Task

**Files:**
- Modify: `mecfs_bio/build_system/task/genome_reference_harmonization/fasta.py`
- Create: `mecfs_bio/build_system/task/genome_reference_harmonization/panel_batch_checks.py`
- Create: `mecfs_bio/build_system/task/genome_reference_harmonization/pan_ukbb/__init__.py` (empty)
- Create: `mecfs_bio/build_system/task/genome_reference_harmonization/pan_ukbb/pan_ukbb_allele_frequency_panel_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_pan_ukbb_allele_frequency_panel_task.py`

**Interfaces:**
- Consumes: `panel_af_col`, `PanUkbbGroup`, `PanelAlleleFrequencyColumns.prefixed` (Task 1); `parquet_encoding`, `open_parquet_writer` (Task 2); `load_fasta` from genome_reference_harmonization_task; `write_fasta` fixture (Task 1).
- Produces (used by Tasks 4 and 5):
  - fasta.py: `gwaslab_code_to_contig_name(code: int) -> str`; `reference_is_acgt(fasta, chrom, positions: np.ndarray, lengths: np.ndarray, max_gather_bytes=DEFAULT_MAX_GATHER_BYTES) -> np.ndarray`.
  - panel_batch_checks.py: `SOURCE_CONTIG_COL = "contig"`; `PanelBatchCounts(rows_in, fasta_ambiguous, ref_mismatch)`; `PanelBatchContext(fasta: IndexedFasta, contig_codes: Mapping[str, int], max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES)`; `PanelWriteSummary(rows_written: int, counts: PanelBatchCounts, chromosomes: frozenset[int])`; `write_checked_panel(batches: Iterable[pl.DataFrame], out_path: Path, context: PanelBatchContext, byte_stream_split_columns: Sequence[str]) -> PanelWriteSummary`; `fetch_file_path(fetch: Fetch, task: Task) -> Path`; `DEFAULT_BATCH_ROWS = 1_000_000`.
  - Panel sources are read with polars (`scan_csv(...).collect_batches(chunk_size=batch_rows)`), which detects compression itself; no reader assumes a compression format.
  - Input batches carry SOURCE_CONTIG_COL (String), POS, REF, ALT and value columns; output tables carry CHR (Int32) first, then POS (Int32), REF, ALT, then the value columns in input order.
  - `PanUkbbAlleleFrequencyPanelTask.create(asset_id, manifest_task, fasta_task, groups, chromosomes, expected_ref_mismatches)`; `PAN_UKBB_MANIFEST_READ_SPEC` (the manifest's DataFrameReadSpec, used by the download asset in Task 6). The manifest is read through `scan_dataframe_asset`, the repo's standard reader, so the download meta must carry that read_spec.

- [ ] **Step 1: Write the failing Pan-UKBB tests**

`test_pan_ukbb_allele_frequency_panel_task.py`:

```python
"""Task-level tests of PanUkbbAlleleFrequencyPanelTask on a synthetic bgzipped manifest."""

from pathlib import Path, PurePath

import attrs
import polars as pl
import pysam
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    PAN_UKBB_MANIFEST_READ_SPEC,
    PanUkbbAlleleFrequencyPanelTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanUkbbGroup,
    panel_af_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    write_fasta,
)

# The manifest's real header (full_variant_qc_metrics.txt.bgz, 2020-08-28).
_HEADER = (
    "chrom pos ref alt rsid varid pass_gnomad_genomes n_passing_populations high_quality "
    "nearest_genes info ac_AFR af_AFR an_AFR gnomad_genomes_ac_AFR gnomad_genomes_af_AFR "
    "gnomad_genomes_an_AFR ac_AMR af_AMR an_AMR gnomad_genomes_ac_AMR gnomad_genomes_af_AMR "
    "gnomad_genomes_an_AMR ac_CSA af_CSA an_CSA ac_EAS af_EAS an_EAS gnomad_genomes_ac_EAS "
    "gnomad_genomes_af_EAS gnomad_genomes_an_EAS ac_EUR af_EUR an_EUR gnomad_genomes_ac_EUR "
    "gnomad_genomes_af_EUR gnomad_genomes_an_EUR ac_MID af_MID an_MID"
).split()
_GROUPS: tuple[PanUkbbGroup, ...] = (
    "ukb_afr",
    "ukb_amr",
    "ukb_csa",
    "ukb_eas",
    "ukb_eur",
    "ukb_mid",
)
# chr1 cycles ACGT and ends in ten Ns (1991-2000); chrX is 200 bases.
_CHR1 = ("ACGT" * 500)[:1990] + "N" * 10
_CHRX = "TTGCA" * 40
_NEXT_BASE = {"A": "C", "C": "G", "G": "T", "T": "A"}
_SWAPPED_POSITIONS = (10, 20)
_BGZF_BLOCK_BYTES = 65536
_MANIFEST_ID = "manifest"
_FASTA_ID = "fasta"


@frozen(slots=True)
class ManifestRow:
    chrom: str
    pos: int
    ref: str
    alt: str
    af: float


def _reference_rows(chrom: str, sequence: str, n_positions: int) -> list[ManifestRow]:
    return [
        ManifestRow(
            chrom=chrom,
            pos=pos,
            ref=sequence[pos - 1],
            alt=_NEXT_BASE[sequence[pos - 1]],
            af=pos / 10_000,
        )
        for pos in range(1, n_positions + 1)
    ]


def _swapped(row: ManifestRow) -> ManifestRow:
    return attrs.evolve(row, ref=row.alt, alt=row.ref)


def _chr1_rows() -> list[ManifestRow]:
    return [
        _swapped(row) if row.pos in _SWAPPED_POSITIONS else row
        for row in _reference_rows("1", _CHR1, 1500)
    ]


def _chrx_rows() -> list[ManifestRow]:
    return _reference_rows("X", _CHRX, 100)


def _write_manifest(path: Path, rows: list[ManifestRow], header: list[str]) -> None:
    lines = ["\t".join(header)]
    for row in rows:
        values = {column: "NA" for column in header} | {
            "chrom": row.chrom,
            "pos": str(row.pos),
            "ref": row.ref,
            "alt": row.alt,
            "varid": f"{row.chrom}:{row.pos}_{row.ref}_{row.alt}",
            "info": "9.5000e-01",
            "nearest_genes": "GENE1,GENE2",
        }
        for column in header:
            if column.startswith("af_"):
                values[column] = f"{row.af:.4e}"
        lines.append("\t".join(values[column] for column in header))
    plain = path.parent / (path.name + ".plain")
    plain.write_text("\n".join(lines) + "\n")
    # Several BGZF blocks, so a reader that stops after the first gzip member is caught.
    assert plain.stat().st_size > 3 * _BGZF_BLOCK_BYTES
    pysam.tabix_compress(str(plain), str(path), force=True)


# Small batches, so the manifest arrives in many batches and chromosome changes and
# swapped rows fall inside and across batch boundaries.
_TEST_BATCH_ROWS = 97


def _run(
    tmp_path: Path,
    rows: list[ManifestRow],
    expected_ref_mismatches: int = len(_SWAPPED_POSITIONS),
    header: list[str] = _HEADER,
) -> pl.DataFrame:
    manifest_path = tmp_path / "manifest.txt.bgz"
    _write_manifest(manifest_path, rows, header)
    fasta_dir = tmp_path / "fasta"
    write_fasta(fasta_dir, {"chr1": _CHR1, "chrX": _CHRX})
    created = PanUkbbAlleleFrequencyPanelTask.create(
        asset_id="pan_ukbb_panel",
        manifest_task=FakeTask(
            ReferenceFileMeta(
                group="pan_ukbb",
                sub_group="variant_manifest",
                sub_folder=PurePath("raw"),
                extension=".txt.bgz",
                id=AssetId(_MANIFEST_ID),
                read_spec=PAN_UKBB_MANIFEST_READ_SPEC,
            )
        ),
        fasta_task=FakeTask(
            FASTAMeta(
                group="genome_sequence",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                id=AssetId(_FASTA_ID),
                build="19",
            )
        ),
        groups=_GROUPS,
        chromosomes=(1, 23),
        expected_ref_mismatches=expected_ref_mismatches,
    )
    task = attrs.evolve(created, batch_rows=_TEST_BATCH_ROWS)
    assets: dict[str, Asset] = {
        _MANIFEST_ID: FileAsset(manifest_path),
        _FASTA_ID: DirectoryAsset(fasta_dir),
    }

    def fetch(asset_id: AssetId) -> Asset:
        return assets[asset_id]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_manifest_becomes_a_panel_without_the_swapped_rows(tmp_path: Path) -> None:
    rows = _chr1_rows() + _chrx_rows()
    panel = _run(tmp_path, rows)
    assert panel.height == len(rows) - len(_SWAPPED_POSITIONS)
    assert panel.columns == [
        GWASLAB_CHROM_COL,
        GWASLAB_POS_COL,
        PANEL_REF_COL,
        PANEL_ALT_COL,
        *[panel_af_col(group) for group in _GROUPS],
    ]
    assert panel.schema[panel_af_col("ukb_eur")] == pl.Float32
    assert panel[GWASLAB_CHROM_COL].unique(maintain_order=True).to_list() == [1, 23]
    chr1 = panel.filter(pl.col(GWASLAB_CHROM_COL) == 1)
    assert not set(_SWAPPED_POSITIONS) & set(chr1[GWASLAB_POS_COL].to_list())
    at_100 = chr1.filter(pl.col(GWASLAB_POS_COL) == 100).row(0, named=True)
    assert at_100[panel_af_col("ukb_eur")] == pytest.approx(0.01)


def test_unexpected_number_of_ref_mismatches_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(
            tmp_path,
            _chr1_rows() + _chrx_rows(),
            expected_ref_mismatches=len(_SWAPPED_POSITIONS) - 1,
        )


def test_row_over_an_n_base_fails(tmp_path: Path) -> None:
    over_n = ManifestRow(chrom="1", pos=1995, ref="A", alt="C", af=0.1)
    with pytest.raises(AssertionError):
        _run(tmp_path, _chr1_rows() + [over_n] + _chrx_rows())


def test_ref_span_beyond_contig_end_fails(tmp_path: Path) -> None:
    beyond = ManifestRow(chrom="X", pos=len(_CHRX) + 1, ref="A", alt="C", af=0.1)
    with pytest.raises(AssertionError):
        _run(tmp_path, _chr1_rows() + _chrx_rows() + [beyond])


def test_missing_frequency_column_fails(tmp_path: Path) -> None:
    header = [column for column in _HEADER if column != "af_MID"]
    with pytest.raises(AssertionError):
        _run(tmp_path, _chr1_rows() + _chrx_rows(), header=header)


def test_chromosomes_out_of_order_fail(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(tmp_path, _chrx_rows() + _chr1_rows())
```

- [ ] **Step 2: Run them to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_pan_ukbb_allele_frequency_panel_task.py -v`
Expected: collection error, ModuleNotFoundError for the pan_ukbb module.

- [ ] **Step 3: Extend fasta.py**

1. Change the constants import to `from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME, GWASLAB_CHROM_NAME_FOR_CODE`, and add `_ACGT_CODES = np.frombuffer(b"ACGT", dtype=np.uint8)` beside `_ASCII_UPPERCASE_MASK`.
2. Add after `contig_to_gwaslab_code`:

```python
def gwaslab_code_to_contig_name(code: int) -> str:
    """Bare contig name (1-22, X, Y, MT) for a gwaslab numeric chromosome code."""
    return GWASLAB_CHROM_NAME_FOR_CODE.get(code, str(code))
```

3. Replace `reference_matches` and `_match_equal_length` with a shared gather, and add `reference_is_acgt`:

```python
@frozen(slots=True)
class _LengthSlice:
    """Row indices whose spans share one length, few enough for one gather."""

    rows: np.ndarray
    length: int


@frozen(slots=True)
class _ReferenceBases:
    """Uppercased reference bases of equal-length spans, and which spans are in bounds."""

    in_bounds: np.ndarray
    bases: np.ndarray

    def __attrs_post_init__(self):
        assert self.in_bounds.ndim == 1 and self.in_bounds.dtype == np.bool_
        assert self.bases.ndim == 2 and self.bases.dtype == np.uint8
        assert self.bases.shape[0] == len(self.in_bounds)


def _length_slices(lengths: np.ndarray, max_gather_bytes: int) -> list[_LengthSlice]:
    """Group rows by span length, split so no gather's offset matrix exceeds the budget."""
    slices: list[_LengthSlice] = []
    order = np.argsort(lengths, kind="stable")
    boundaries = np.flatnonzero(np.diff(lengths[order])) + 1
    for group in np.split(order, boundaries):
        length = int(lengths[group[0]])
        rows_per_slice = max(1, max_gather_bytes // (length * _OFFSET_BYTES))
        for start in range(0, len(group), rows_per_slice):
            slices.append(
                _LengthSlice(rows=group[start : start + rows_per_slice], length=length)
            )
    return slices


def _gather_reference(
    genome: np.ndarray, entry: FaiEntry, positions: np.ndarray, length: int
) -> _ReferenceBases:
    start = positions.astype(np.int64) - 1
    in_bounds = (start >= 0) & (start + length <= entry.length)
    clipped = np.clip(start, 0, max(entry.length - length, 0))
    base_index = clipped[:, None] + np.arange(length, dtype=np.int64)[None, :]
    file_offsets = (
        entry.offset
        + base_index
        + (base_index // entry.line_bases) * (entry.line_bytes - entry.line_bases)
    )
    return _ReferenceBases(
        in_bounds=in_bounds, bases=genome[file_offsets] & _ASCII_UPPERCASE_MASK
    )


def reference_matches(
    fasta: IndexedFasta,
    chrom: int,
    positions: np.ndarray,
    alleles: pl.Series,
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES,
) -> np.ndarray:
    """Boolean array: allele i equals the reference sequence starting at 1-based positions[i]."""
    assert positions.ndim == 1 and len(positions) == len(alleles), (
        "positions and alleles must be one-dimensional and the same length"
    )
    assert chrom in fasta.entries, f"chromosome {chrom} is not in {fasta.fasta_path}"
    assert max_gather_bytes > 0
    result = np.zeros(len(positions), dtype=bool)
    if len(positions) == 0:
        return result
    entry = fasta.entries[chrom]
    genome = np.memmap(fasta.fasta_path, dtype=np.uint8, mode="r")
    lengths = alleles.str.len_bytes().to_numpy()
    assert (lengths > 0).all(), "alleles must be non-empty"
    for piece in _length_slices(lengths, max_gather_bytes):
        reference = _gather_reference(genome, entry, positions[piece.rows], piece.length)
        observed = (
            np.frombuffer(
                "".join(alleles.gather(piece.rows).to_list()).encode("ascii"),
                dtype=np.uint8,
            ).reshape(len(piece.rows), piece.length)
            & _ASCII_UPPERCASE_MASK
        )
        result[piece.rows] = reference.in_bounds & (reference.bases == observed).all(
            axis=1
        )
    return result


def reference_is_acgt(
    fasta: IndexedFasta,
    chrom: int,
    positions: np.ndarray,
    lengths: np.ndarray,
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES,
) -> np.ndarray:
    """Boolean array: the lengths[i] reference bases from 1-based positions[i] are in bounds
    and each is A, C, G or T in either case. False where the span holds N or another IUPAC
    ambiguity code."""
    assert positions.ndim == 1 and lengths.ndim == 1 and len(positions) == len(lengths), (
        "positions and lengths must be one-dimensional and the same length"
    )
    assert chrom in fasta.entries, f"chromosome {chrom} is not in {fasta.fasta_path}"
    assert max_gather_bytes > 0
    result = np.zeros(len(positions), dtype=bool)
    if len(positions) == 0:
        return result
    assert (lengths > 0).all(), "spans must be non-empty"
    entry = fasta.entries[chrom]
    genome = np.memmap(fasta.fasta_path, dtype=np.uint8, mode="r")
    for piece in _length_slices(lengths, max_gather_bytes):
        reference = _gather_reference(genome, entry, positions[piece.rows], piece.length)
        result[piece.rows] = reference.in_bounds & np.isin(
            reference.bases, _ACGT_CODES
        ).all(axis=1)
    return result
```

- [ ] **Step 4: Create panel_batch_checks.py**

```python
"""
Invariants shared by the allele-frequency panel extraction Tasks, applied batch by batch.

A panel source (a gnomAD sites VCF, the Pan-UKBB variant manifest) is streamed as polars
batches holding a source contig column, POS, REF, ALT and frequency columns. Each batch is
checked and converted here, and the surviving rows are written to one parquet file.

Fatal checks:
- contig names must be among the expected ones; they become an Int32 gwaslab CHR code;
- REF and ALT must be non-null and contain only A, C, G and T;
- (CHR, POS) must never decrease, and no (CHR, POS, REF, ALT) key may repeat. The keys at
  the last position of each batch are carried into the next, so a pair split across a batch
  boundary is still caught;
- every REF span must lie inside its FASTA contig (a span past the end means a wrong build or
  contig naming).

REF is then compared with the FASTA. Records whose reference span holds N or another IUPAC
code are dropped and counted as fasta_ambiguous: their orientation is unverifiable, and the
harmonizer's own FASTA classification would reject them anyway. Records whose REF differs
from a pure A/C/G/T span are dropped and counted as ref_mismatch. The calling Task decides
how many of each it accepts.
"""

from collections.abc import Iterable, Mapping, Sequence
from pathlib import Path

import numpy as np
import polars as pl
import pyarrow.parquet
from attrs import frozen

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    open_parquet_writer,
    parquet_encoding,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    DEFAULT_MAX_GATHER_BYTES,
    IndexedFasta,
    reference_is_acgt,
    reference_matches,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

SOURCE_CONTIG_COL = "contig"
PANEL_KEY_COLUMNS = [GWASLAB_CHROM_COL, GWASLAB_POS_COL, PANEL_REF_COL, PANEL_ALT_COL]
# Rows per streamed batch; each batch becomes one parquet row group. On the full Pan-UKBB
# manifest, 1,000,000-row batches peaked at about 3 GB RSS.
DEFAULT_BATCH_ROWS = 1_000_000
_ACGT_ALLELE = r"^[ACGT]+$"
_FASTA_ACGT_COL = "_fasta_acgt"
_REF_MATCHES_COL = "_ref_matches"
_CHROM_STEP_COL = "_chrom_step"
_POS_STEP_COL = "_pos_step"


@frozen(slots=True)
class PanelBatchCounts:
    rows_in: int
    fasta_ambiguous: int
    ref_mismatch: int

    @classmethod
    def zero(cls) -> "PanelBatchCounts":
        return cls(rows_in=0, fasta_ambiguous=0, ref_mismatch=0)

    def __add__(self, other: "PanelBatchCounts") -> "PanelBatchCounts":
        return PanelBatchCounts(
            rows_in=self.rows_in + other.rows_in,
            fasta_ambiguous=self.fasta_ambiguous + other.fasta_ambiguous,
            ref_mismatch=self.ref_mismatch + other.ref_mismatch,
        )


@frozen(slots=True)
class PanelBatchContext:
    """What every batch is checked against: the FASTA and the allowed source contigs."""

    fasta: IndexedFasta
    contig_codes: Mapping[str, int]
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES


@frozen(slots=True)
class CheckedPanelBatch:
    table: pl.DataFrame
    counts: PanelBatchCounts
    carry: pl.DataFrame


@frozen(slots=True)
class PanelWriteSummary:
    rows_written: int
    counts: PanelBatchCounts
    chromosomes: frozenset[int]


def fetch_file_path(fetch: Fetch, task: Task) -> Path:
    """Fetch a Task's FileAsset and return its path."""
    asset = fetch(task.asset_id)
    assert isinstance(asset, FileAsset), f"expected {task.asset_id} to be a FileAsset"
    return asset.path


def write_checked_panel(
    batches: Iterable[pl.DataFrame],
    out_path: Path,
    context: PanelBatchContext,
    byte_stream_split_columns: Sequence[str],
) -> PanelWriteSummary:
    """Check every batch, write the kept rows to out_path as one parquet, and summarize."""
    counts = PanelBatchCounts.zero()
    carry = _empty_keys()
    rows_written = 0
    chromosomes: set[int] = set()
    writer: pyarrow.parquet.ParquetWriter | None = None
    try:
        for batch in batches:
            checked = check_panel_batch(batch, context, carry)
            carry = checked.carry
            counts = counts + checked.counts
            table = checked.table.to_arrow()
            if writer is None:
                writer = open_parquet_writer(
                    out_path,
                    table.schema,
                    parquet_encoding(
                        table.schema.names,
                        compression="zstd",
                        compression_level=None,
                        byte_stream_split_columns=byte_stream_split_columns,
                    ),
                )
            writer.write_table(table)
            rows_written += checked.table.height
            chromosomes.update(checked.table[GWASLAB_CHROM_COL].unique().to_list())
    finally:
        if writer is not None:
            writer.close()
    assert writer is not None, "the panel source produced no batches"
    return PanelWriteSummary(
        rows_written=rows_written, counts=counts, chromosomes=frozenset(chromosomes)
    )


def check_panel_batch(
    batch: pl.DataFrame, context: PanelBatchContext, carry: pl.DataFrame
) -> CheckedPanelBatch:
    """Apply the module's checks to one batch; carry holds the previous batch's last keys."""
    converted = _with_chromosome_codes(batch, context.contig_codes)
    _assert_acgt_alleles(converted)
    keys = pl.concat([carry, converted.select(PANEL_KEY_COLUMNS)])
    _assert_sorted_and_unique(keys)
    classified = _classify_against_fasta(converted, context)
    ambiguous = classified.select((~pl.col(_FASTA_ACGT_COL)).sum()).item()
    mismatched = classified.select(
        (pl.col(_FASTA_ACGT_COL) & ~pl.col(_REF_MATCHES_COL)).sum()
    ).item()
    return CheckedPanelBatch(
        table=classified.filter(
            pl.col(_FASTA_ACGT_COL) & pl.col(_REF_MATCHES_COL)
        ).drop(_FASTA_ACGT_COL, _REF_MATCHES_COL),
        counts=PanelBatchCounts(
            rows_in=converted.height,
            fasta_ambiguous=int(ambiguous),
            ref_mismatch=int(mismatched),
        ),
        carry=_keys_at_last_position(keys, carry),
    )


def _empty_keys() -> pl.DataFrame:
    return pl.DataFrame(
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
        }
    )


def _with_chromosome_codes(
    batch: pl.DataFrame, contig_codes: Mapping[str, int]
) -> pl.DataFrame:
    contigs = set(batch[SOURCE_CONTIG_COL].unique().to_list())
    unexpected = contigs - set(contig_codes)
    assert not unexpected, (
        f"unexpected contigs {sorted(map(str, unexpected))}; expected {sorted(contig_codes)}"
    )
    value_columns = [
        column
        for column in batch.columns
        if column
        not in (SOURCE_CONTIG_COL, GWASLAB_POS_COL, PANEL_REF_COL, PANEL_ALT_COL)
    ]
    return batch.select(
        pl.col(SOURCE_CONTIG_COL)
        .replace_strict(dict(contig_codes), return_dtype=pl.Int32)
        .alias(GWASLAB_CHROM_COL),
        pl.col(GWASLAB_POS_COL).cast(pl.Int32),
        PANEL_REF_COL,
        PANEL_ALT_COL,
        *value_columns,
    )


def _assert_acgt_alleles(batch: pl.DataFrame) -> None:
    for column in (PANEL_REF_COL, PANEL_ALT_COL):
        bad = batch.filter(
            pl.col(column).is_null() | ~pl.col(column).str.contains(_ACGT_ALLELE)
        )
        assert bad.height == 0, (
            f"{bad.height} records with a null or non-ACGT {column}, e.g. "
            f"{bad.select(PANEL_KEY_COLUMNS).head(3).rows()}"
        )


def _assert_sorted_and_unique(keys: pl.DataFrame) -> None:
    steps = keys.with_columns(
        pl.col(GWASLAB_CHROM_COL).diff().alias(_CHROM_STEP_COL),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64).diff().alias(_POS_STEP_COL),
    )
    out_of_order = steps.filter(
        (pl.col(_CHROM_STEP_COL) < 0)
        | ((pl.col(_CHROM_STEP_COL) == 0) & (pl.col(_POS_STEP_COL) < 0))
    )
    assert out_of_order.height == 0, (
        "records out of (CHR, POS) order at "
        f"{out_of_order.select(GWASLAB_CHROM_COL, GWASLAB_POS_COL).head(3).rows()}"
    )
    duplicated = keys.filter(pl.struct(PANEL_KEY_COLUMNS).is_duplicated())
    assert duplicated.height == 0, (
        f"duplicate panel keys, e.g. {duplicated.unique().head(3).rows()}"
    )


def _keys_at_last_position(keys: pl.DataFrame, carry: pl.DataFrame) -> pl.DataFrame:
    if keys.height == 0:
        return carry
    last = keys.row(-1, named=True)
    return keys.filter(
        (pl.col(GWASLAB_CHROM_COL) == last[GWASLAB_CHROM_COL])
        & (pl.col(GWASLAB_POS_COL) == last[GWASLAB_POS_COL])
    )


def _classify_against_fasta(
    batch: pl.DataFrame, context: PanelBatchContext
) -> pl.DataFrame:
    if batch.height == 0:
        return batch.with_columns(
            pl.lit(True).alias(_FASTA_ACGT_COL), pl.lit(True).alias(_REF_MATCHES_COL)
        )
    return pl.concat(
        [
            _classify_chromosome(rows, int(code), context)
            for (code,), rows in batch.partition_by(
                GWASLAB_CHROM_COL, as_dict=True, maintain_order=True
            ).items()
        ]
    )


def _classify_chromosome(
    rows: pl.DataFrame, code: int, context: PanelBatchContext
) -> pl.DataFrame:
    fasta = context.fasta
    assert code in fasta.entries, f"chromosome {code} is not in {fasta.fasta_path}"
    positions = rows[GWASLAB_POS_COL].to_numpy()
    lengths = rows[PANEL_REF_COL].str.len_bytes().to_numpy()
    span_end = positions.astype(np.int64) + lengths - 1
    beyond = (positions < 1) | (span_end > fasta.entries[code].length)
    assert not beyond.any(), (
        f"chromosome {code}: {int(beyond.sum())} REF spans lie outside the FASTA contig "
        "(wrong build or contig naming?)"
    )
    return rows.with_columns(
        pl.Series(
            _FASTA_ACGT_COL,
            reference_is_acgt(
                fasta, code, positions, lengths, context.max_gather_bytes
            ),
        ),
        pl.Series(
            _REF_MATCHES_COL,
            reference_matches(
                fasta, code, positions, rows[PANEL_REF_COL], context.max_gather_bytes
            ),
        ),
    )
```

- [ ] **Step 5: Create the Pan-UKBB Task**

`pan_ukbb/pan_ukbb_allele_frequency_panel_task.py`:

```python
"""
Allele-frequency panel from the Pan-UK Biobank variant manifest, for genome-reference
harmonization.

The manifest (full_variant_qc_metrics.txt.bgz) is a bgzipped TSV of UK Biobank imputed
variants on GRCh37, already filtered to imputation INFO > 0.8, with the alternate allele
frequency af_{POP} for the groups AFR, AMR, CSA, EAS, EUR and MID. This Task streams it once
and writes CHR, POS, REF, ALT and one AF_ukb_{pop} column per configured group, through the
shared panel batch checks.

Every manifest row is kept except those whose ref disagrees with the FASTA. The manifest
holds a known set of such rows: ref and alt swapped, inherited from UK Biobank's imputed
BGEN allele order. expected_ref_mismatches pins their number, so any change in the file or in
the parsing fails the build. Allele numbers are not stored: in the manifest they are
constant within a contig and carry no per-site information.

The manifest is read through the asset's read_spec with polars, which detects the
compression, and streamed in batches of batch_rows rows. collect_batches is marked unstable
in polars; the Task tests pin that every row of a multi-block bgzipped file arrives.

See experiments/claude/design_specs/2026-10-02-gnomad-allele-frequency-panel-design.md.
"""

from collections.abc import Iterator, Sequence
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    gwaslab_code_to_contig_name,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    load_fasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    DEFAULT_BATCH_ROWS,
    SOURCE_CONTIG_COL,
    PanelBatchContext,
    PanelWriteSummary,
    write_checked_panel,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanUkbbGroup,
    panel_af_col,
)
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL

logger = structlog.get_logger()

PAN_UKBB_BUILD: GenomeBuild = "19"
PAN_UKBB_PANEL_FILENAME = "pan_ukbb_allele_frequencies"
_MANIFEST_CONTIG_COL = "chrom"
_MANIFEST_POS_COL = "pos"
_MANIFEST_REF_COL = "ref"
_MANIFEST_ALT_COL = "alt"
# How to read the manifest. chrom must be read as a string: schema inference over the
# leading rows would otherwise type it as an integer and fail at "X".
PAN_UKBB_MANIFEST_READ_SPEC = DataFrameReadSpec(
    DataFrameTextFormat(
        separator="\t",
        null_values=["NA"],
        schema_overrides={_MANIFEST_CONTIG_COL: pl.String()},
    )
)


def manifest_af_col(group: PanUkbbGroup) -> str:
    """The manifest's frequency column for a group: ukb_eur -> af_EUR."""
    return "af_" + group.removeprefix("ukb_").upper()


@frozen(slots=True)
class PanUkbbAlleleFrequencyPanelTask(Task):
    meta: HarmonizableReferenceTableMeta
    manifest_task: Task
    fasta_task: Task
    groups: tuple[PanUkbbGroup, ...]
    chromosomes: tuple[int, ...]
    expected_ref_mismatches: int
    batch_rows: int = DEFAULT_BATCH_ROWS

    def __attrs_post_init__(self):
        assert self.groups, "at least one Pan-UKBB group is required"
        assert len(set(self.groups)) == len(self.groups), f"duplicate groups {self.groups}"
        assert self.chromosomes, "at least one chromosome is required"
        assert self.expected_ref_mismatches >= 0
        assert self.batch_rows > 0

    @property
    def deps(self) -> list[Task]:
        return [self.manifest_task, self.fasta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        manifest = scan_manifest(
            fetch(self.manifest_task.asset_id), self.manifest_task.meta
        )
        fasta = load_fasta(fetch, self.fasta_task)
        assert_manifest_columns(manifest, self.groups)
        out_path = scratch_dir / (PAN_UKBB_PANEL_FILENAME + ".parquet")
        summary = write_checked_panel(
            manifest_batches(manifest, self.groups, self.batch_rows),
            out_path,
            PanelBatchContext(
                fasta=fasta, contig_codes=manifest_contig_codes(self.chromosomes)
            ),
            byte_stream_split_columns=[],
        )
        assert_pan_ukbb_summary(
            summary, self.expected_ref_mismatches, self.chromosomes
        )
        logger.info(
            "Pan-UKBB allele-frequency panel written",
            rows_in=summary.counts.rows_in,
            rows_written=summary.rows_written,
            ref_mismatches_dropped=summary.counts.ref_mismatch,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        manifest_task: Task,
        fasta_task: Task,
        groups: Sequence[PanUkbbGroup],
        chromosomes: Sequence[int],
        expected_ref_mismatches: int,
    ) -> "PanUkbbAlleleFrequencyPanelTask":
        fasta_meta = fasta_task.meta
        assert isinstance(fasta_meta, FASTAMeta), (
            f"fasta_task must carry FASTAMeta, got {type(fasta_meta).__name__}"
        )
        assert fasta_meta.build == PAN_UKBB_BUILD, (
            f"the Pan-UKBB manifest is GRCh37; got a build {fasta_meta.build} FASTA"
        )
        source_meta = manifest_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta manifest, got {type(source_meta).__name__}"
        )
        assert source_meta.read_spec is not None, (
            "the manifest meta needs a read_spec (PAN_UKBB_MANIFEST_READ_SPEC)"
        )
        return cls(
            meta=HarmonizableReferenceTableMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename=PAN_UKBB_PANEL_FILENAME,
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=PAN_UKBB_BUILD,
                    ref_allele_col=PANEL_REF_COL,
                    pos_col=GWASLAB_POS_COL,
                ),
                allele_frequency_columns=PanelAlleleFrequencyColumns.prefixed(groups),
            ),
            manifest_task=manifest_task,
            fasta_task=fasta_task,
            groups=tuple(groups),
            chromosomes=tuple(chromosomes),
            expected_ref_mismatches=expected_ref_mismatches,
        )


def manifest_contig_codes(chromosomes: Sequence[int]) -> dict[str, int]:
    """Manifest contig names ("1".."22", "X") for the expected gwaslab codes."""
    return {gwaslab_code_to_contig_name(code): code for code in chromosomes}


def scan_manifest(asset: Asset, meta: Meta) -> pl.LazyFrame:
    """The manifest as a polars LazyFrame, read through the asset's read_spec."""
    native = scan_dataframe_asset(asset, meta).to_native()
    assert isinstance(native, pl.LazyFrame), (
        f"expected a polars LazyFrame for the manifest, got {type(native).__name__}"
    )
    return native


def assert_manifest_columns(
    manifest: pl.LazyFrame, groups: Sequence[PanUkbbGroup]
) -> None:
    present = manifest.collect_schema().names()
    required = [
        _MANIFEST_CONTIG_COL,
        _MANIFEST_POS_COL,
        _MANIFEST_REF_COL,
        _MANIFEST_ALT_COL,
        *[manifest_af_col(group) for group in groups],
    ]
    missing = [column for column in required if column not in present]
    assert not missing, f"the Pan-UKBB manifest lacks columns {missing}"


def manifest_batches(
    manifest: pl.LazyFrame, groups: Sequence[PanUkbbGroup], batch_rows: int
) -> Iterator[pl.DataFrame]:
    """Stream the manifest's needed columns, renamed to the panel's, in batches."""
    af_columns = {manifest_af_col(group): panel_af_col(group) for group in groups}
    query = manifest.select(
        pl.col(_MANIFEST_CONTIG_COL).cast(pl.String).alias(SOURCE_CONTIG_COL),
        pl.col(_MANIFEST_POS_COL).alias(GWASLAB_POS_COL),
        pl.col(_MANIFEST_REF_COL).alias(PANEL_REF_COL),
        pl.col(_MANIFEST_ALT_COL).alias(PANEL_ALT_COL),
        *[
            pl.col(source).cast(pl.Float32).alias(target)
            for source, target in af_columns.items()
        ],
    )
    for batch in query.collect_batches(chunk_size=batch_rows):
        null_counts = batch.select(
            [pl.col(column).null_count() for column in af_columns.values()]
        ).row(0, named=True)
        assert not any(null_counts.values()), (
            f"null allele frequencies in the Pan-UKBB manifest: {null_counts}"
        )
        yield batch


def assert_pan_ukbb_summary(
    summary: PanelWriteSummary,
    expected_ref_mismatches: int,
    chromosomes: Sequence[int],
) -> None:
    assert summary.counts.ref_mismatch == expected_ref_mismatches, (
        f"expected {expected_ref_mismatches} manifest rows whose ref differs from the "
        f"FASTA, found {summary.counts.ref_mismatch}"
    )
    assert summary.counts.fasta_ambiguous == 0, (
        f"{summary.counts.fasta_ambiguous} manifest rows lie over N or IUPAC FASTA bases"
    )
    missing = sorted(set(chromosomes) - summary.chromosomes)
    assert not missing, f"no Pan-UKBB rows on chromosomes {missing}"
    assert summary.rows_written > 0, "the Pan-UKBB panel is empty"
```

- [ ] **Step 6: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/ -v`
Expected: all PASS (the six new tests, and the harmonization tests, which exercise the refactored reference_matches).

- [ ] **Step 7: Run green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_task3.log
git add mecfs_bio/build_system/task/genome_reference_harmonization test_mecfs_bio/unit/build_system/task/genome_reference_harmonization
git commit -m "Add shared panel batch checks and the Pan-UKBB panel Task

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: gnomAD release descriptor and per-chromosome Task

**Files:**
- Create: `mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/__init__.py` (empty)
- Create: `.../gnomad/gnomad_release.py`
- Create: `.../gnomad/gnomad_chromosome_allele_frequency_task.py`
- Test: `test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_chromosome_allele_frequency_task.py`

**Interfaces:**
- Consumes: Task 3's panel_batch_checks API, `gwaslab_code_to_contig_name`; Task 1's `panel_af_col`, `panel_an_col`, `GnomadGroup`; `load_fasta`; `execute_command_with_retries` from `mecfs_bio.util.subproc.run_command`.
- Produces (used by Tasks 5 and 6):
  - `GnomadRelease(name, build, vcf_url_template, contig_prefix: ContigPrefix, chromosomes, main_groups, extra_groups, header_assembly)` with property `groups`; `CHROM_PLACEHOLDER = "{chrom}"`; `gnomad_contig_name(release, chrom) -> str`; `gnomad_vcf_url(release, chrom) -> str`; `GNOMAD_GROUP = "gnomad"`.
  - `GnomadChromosomeAlleleFrequencyTask(meta: ReferenceFileMeta, release, chrom, fasta_task, batch_rows=DEFAULT_BATCH_ROWS, max_attempts=DEFAULT_MAX_ATTEMPTS)`; `.create(release, chrom, fasta_task)`; `read_vcf_header(url, max_attempts) -> str`; `assert_gnomad_header(header, release, chrom)`.
  - Output parquet columns: CHR, POS, REF, ALT, then AF_g and AN_g for each group of release.groups, interleaved (AF_g1, AN_g1, AF_g2, AN_g2, ...).

- [ ] **Step 1: Write the failing tests**

`test_gnomad_chromosome_allele_frequency_task.py`:

```python
"""Task-level tests of GnomadChromosomeAlleleFrequencyTask on tiny local bgzipped VCFs."""

from pathlib import Path, PurePath
from subprocess import CalledProcessError

import attrs
import polars as pl
import pyarrow.parquet as pq
import pysam
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    DEFAULT_BATCH_ROWS,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    GnomadGroup,
    panel_af_col,
    panel_an_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    write_fasta,
)

# chr1: ACGT repeated (pos p holds "ACGT"[(p - 1) % 4]), with N at 101-104.
_CHR1 = "ACGT" * 25 + "NNNN" + "ACGT" * 24
_CHRX = "TTGCA" * 20
_MAIN_GROUPS: tuple[GnomadGroup, ...] = ("afr", "nfe")
_EXTRA_GROUPS: tuple[GnomadGroup, ...] = ("nfe_nwe",)
_ALL_GROUPS = _MAIN_GROUPS + _EXTRA_GROUPS
_ASSEMBLY = "gnomAD_GRCh37"
_FASTA_ID = "fasta"
# One row per batch, so every adjacent pair of records straddles a batch boundary.
_ONE_ROW_BATCH = 1
_LONG_DELETION_REF = _CHR1[8:38]  # 30 bases from pos 9


@frozen(slots=True)
class Record:
    pos: int
    ref: str
    alt: str
    filter: str = "PASS"
    afr: str = "0"
    afr_an: int = 100
    nfe: str = "0"
    nfe_an: int = 200
    nfe_nwe: str = "0"
    nfe_nwe_an: int = 100


def _header(groups: tuple[GnomadGroup, ...] = _ALL_GROUPS) -> list[str]:
    lines = [
        "##fileformat=VCFv4.2",
        '##FILTER=<ID=PASS,Description="All filters passed">',
        '##FILTER=<ID=RF,Description="Failed random forest">',
        f"##contig=<ID=1,length={len(_CHR1)},assembly={_ASSEMBLY}>",
        f"##contig=<ID=X,length={len(_CHRX)},assembly={_ASSEMBLY}>",
    ]
    for group in groups:
        lines.append(
            f'##INFO=<ID={panel_af_col(group)},Number=A,Type=Float,Description="AF {group}">'
        )
        lines.append(
            f'##INFO=<ID={panel_an_col(group)},Number=1,Type=Integer,Description="AN {group}">'
        )
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")
    return lines


def _line(record: Record) -> str:
    info = ";".join(
        [
            f"AF_afr={record.afr};AN_afr={record.afr_an}",
            f"AF_nfe={record.nfe};AN_nfe={record.nfe_an}",
            f"AF_nfe_nwe={record.nfe_nwe};AN_nfe_nwe={record.nfe_nwe_an}",
        ]
    )
    return f"1\t{record.pos}\t.\t{record.ref}\t{record.alt}\t.\t{record.filter}\t{info}"


def _release(tmp_path: Path) -> GnomadRelease:
    return GnomadRelease(
        name="test_release",
        build="19",
        vcf_url_template=str(tmp_path / "sites.{chrom}.vcf.bgz"),
        contig_prefix="",
        chromosomes=(1, 23),
        main_groups=_MAIN_GROUPS,
        extra_groups=_EXTRA_GROUPS,
        header_assembly=_ASSEMBLY,
    )


def _run(
    tmp_path: Path,
    records: list[Record],
    header: list[str] | None = None,
    batch_rows: int = DEFAULT_BATCH_ROWS,
    truncate: bool = False,
) -> Path:
    plain = tmp_path / "sites.1.vcf"
    plain.write_text(
        "\n".join((header or _header()) + [_line(r) for r in records]) + "\n"
    )
    vcf = tmp_path / "sites.1.vcf.bgz"
    pysam.tabix_compress(str(plain), str(vcf), force=True)
    if truncate:
        vcf.write_bytes(vcf.read_bytes()[:-40])
    fasta_dir = tmp_path / "fasta"
    write_fasta(fasta_dir, {"chr1": _CHR1, "chrX": _CHRX})
    fasta_task = FakeTask(
        FASTAMeta(
            group="genome_sequence",
            sub_group="synthetic",
            sub_folder=PurePath("processed"),
            id=AssetId(_FASTA_ID),
            build="19",
        )
    )
    task = attrs.evolve(
        GnomadChromosomeAlleleFrequencyTask.create(
            release=_release(tmp_path), chrom=1, fasta_task=fasta_task
        ),
        batch_rows=batch_rows,
        max_attempts=1,
    )

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == _FASTA_ID
        return DirectoryAsset(fasta_dir)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return result.path


def _encodings(path: Path, column: str) -> set[str]:
    row_group = pq.ParquetFile(path).metadata.row_group(0)
    columns = [row_group.column(i) for i in range(row_group.num_columns)]
    return set(next(c for c in columns if c.path_in_schema == column).encodings)


_HAPPY_RECORDS = [
    Record(1, "A", "G", afr="0.1"),  # kept
    Record(1, "A", "T", nfe="0.01"),  # same site, other ALT: kept, not a duplicate
    Record(2, "C", "T", filter="RF", afr="0.1"),  # not PASS: dropped
    Record(3, "G", "A", nfe_nwe="0.01"),  # polymorphic only in an extra group: dropped
    Record(4, "T", "C", afr=".", afr_an=0, nfe="0.2"),  # kept, AF_afr null (AN 0)
    Record(5, "A", "AC", nfe="0.3"),  # insertion: kept
    Record(9, _LONG_DELETION_REF, "A", nfe="0.05"),  # 30-base deletion: kept
    Record(102, "C", "G", afr="0.1"),  # REF over N: dropped and counted, not fatal
]


@pytest.mark.parametrize("batch_rows", [DEFAULT_BATCH_ROWS, _ONE_ROW_BATCH])
def test_polymorphic_pass_records_become_an_allele_frequency_table(
    tmp_path: Path, batch_rows: int
) -> None:
    path = _run(tmp_path, _HAPPY_RECORDS, batch_rows=batch_rows)
    table = pl.read_parquet(path)
    assert table.select(GWASLAB_POS_COL, "ALT").rows() == [
        (1, "G"),
        (1, "T"),
        (4, "C"),
        (5, "AC"),
        (9, "A"),
    ]
    assert table[GWASLAB_CHROM_COL].unique().to_list() == [1]
    assert table.schema[GWASLAB_CHROM_COL] == pl.Int32
    assert table.schema[panel_af_col("afr")] == pl.Float32
    assert table.schema[panel_an_col("afr")] == pl.Int32
    at_4 = table.filter(pl.col(GWASLAB_POS_COL) == 4).row(0, named=True)
    assert at_4[panel_af_col("afr")] is None
    assert at_4[panel_an_col("afr")] == 0
    for group in _ALL_GROUPS:
        assert "BYTE_STREAM_SPLIT" in _encodings(path, panel_an_col(group))
        assert "RLE_DICTIONARY" in _encodings(path, panel_af_col(group))


def test_ref_mismatch_over_acgt_reference_fails(tmp_path: Path) -> None:
    # pos 6 holds C
    with pytest.raises(AssertionError):
        _run(tmp_path, [Record(6, "A", "G", afr="0.1")])


def test_header_missing_a_group_field_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(
            tmp_path,
            [Record(1, "A", "G", afr="0.1")],
            header=_header(groups=("afr", "nfe")),
        )


@pytest.mark.parametrize("batch_rows", [DEFAULT_BATCH_ROWS, _ONE_ROW_BATCH])
def test_duplicate_key_fails(tmp_path: Path, batch_rows: int) -> None:
    records = [Record(1, "A", "G", afr="0.1"), Record(1, "A", "G", afr="0.2")]
    with pytest.raises(AssertionError):
        _run(tmp_path, records, batch_rows=batch_rows)


def test_decreasing_position_fails(tmp_path: Path) -> None:
    records = [Record(5, "A", "G", afr="0.1"), Record(1, "A", "G", afr="0.1")]
    with pytest.raises(AssertionError):
        _run(tmp_path, records)


def test_chromosome_without_polymorphic_pass_records_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(tmp_path, [Record(1, "A", "G", filter="RF", afr="0.1")])


def test_truncated_vcf_fails(tmp_path: Path) -> None:
    with pytest.raises(CalledProcessError):
        _run(tmp_path, _HAPPY_RECORDS, truncate=True)
```

- [ ] **Step 2: Run them to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_chromosome_allele_frequency_task.py -v`
Expected: collection error, ModuleNotFoundError for the gnomad package.

- [ ] **Step 3: Create gnomad_release.py**

```python
"""
Descriptor of one gnomAD sites-VCF release: where its per-chromosome files live and which
ancestry groups it carries. Instances live under mecfs_bio/assets/reference_data/gnomad/.
"""

from typing import Literal

from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    gwaslab_code_to_contig_name,
)
from mecfs_bio.constants.allele_frequency_panel_constants import GnomadGroup
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME

ContigPrefix = Literal["", "chr"]
CHROM_PLACEHOLDER = "{chrom}"
# Asset-store group shared by every gnomAD asset.
GNOMAD_GROUP = "gnomad"


@frozen(slots=True)
class GnomadRelease:
    """A gnomAD release's per-chromosome sites VCFs.

    vcf_url_template contains CHROM_PLACEHOLDER, replaced by the bare chromosome name
    (1-22, X, Y). contig_prefix is what the VCF prepends to contig names ("chr" in GRCh38
    releases). main_groups admit records (a record is kept when polymorphic in one of
    them); extra_groups are stored but admit nothing. header_assembly is the assembly
    string every ##contig header line must carry.
    """

    name: str
    build: GenomeBuild
    vcf_url_template: str
    contig_prefix: ContigPrefix
    chromosomes: tuple[int, ...]
    main_groups: tuple[GnomadGroup, ...]
    extra_groups: tuple[GnomadGroup, ...]
    header_assembly: str

    def __attrs_post_init__(self):
        assert CHROM_PLACEHOLDER in self.vcf_url_template, (
            f"vcf_url_template must contain {CHROM_PLACEHOLDER}"
        )
        assert self.main_groups, "a release needs at least one main group"
        assert len(set(self.groups)) == len(self.groups), (
            f"groups must be distinct, and main and extra disjoint: {self.groups}"
        )
        assert self.chromosomes, "a release needs at least one chromosome"
        assert list(self.chromosomes) == sorted(set(self.chromosomes)), (
            f"chromosomes must be ascending gwaslab codes without repeats: {self.chromosomes}"
        )
        assert GWASLAB_CHROM_CODE_FOR_NAME["MT"] not in self.chromosomes, (
            "gnomAD panels exclude MT"
        )

    @property
    def groups(self) -> tuple[GnomadGroup, ...]:
        return self.main_groups + self.extra_groups


def gnomad_contig_name(release: GnomadRelease, chrom: int) -> str:
    """The contig name the release's VCF uses for a gwaslab chromosome code."""
    return release.contig_prefix + gwaslab_code_to_contig_name(chrom)


def gnomad_vcf_url(release: GnomadRelease, chrom: int) -> str:
    return release.vcf_url_template.replace(
        CHROM_PLACEHOLDER, gwaslab_code_to_contig_name(chrom)
    )
```

- [ ] **Step 4: Create gnomad_chromosome_allele_frequency_task.py**

```python
"""
One chromosome of a gnomAD sites VCF, reduced to an allele-frequency table.

Streams the release's per-chromosome VCF (usually over HTTPS) with a single bcftools query
process, keeping PASS records polymorphic in at least one main group, and writes CHR, POS,
REF, ALT and AF_g, AN_g for every group of the release. One Task per chromosome keeps the
genome-wide build resumable: a failure costs one chromosome.

The header is checked first (every AF_g and AN_g declared, every contig line carrying the
release's assembly, this chromosome's contig declared), so a wrong release description fails
in seconds rather than after a long stream. bcftools query runs as one process rather than
a view | query pipe: execute_command runs commands through sh, where a pipeline's exit status
is that of its last command, so an upstream network failure could otherwise leave a
truncated table and exit 0. A single process lets htslib's BGZF integrity checks fail the
Task. Each retry overwrites the scratch TSV.

The TSV becomes parquet in bounded memory through the shared panel batch checks. For every
group, AF must be null exactly where AN is 0. Any REF that differs from a pure A/C/G/T FASTA
span is fatal; records over N or other IUPAC bases are dropped and counted. AN columns are
written with byte-stream-split encoding and AF columns with dictionary encoding, the smallest
combination measured.
"""

import re
import shlex
from collections.abc import Iterator
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    gwaslab_code_to_contig_name,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    load_fasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GNOMAD_GROUP,
    GnomadRelease,
    gnomad_contig_name,
    gnomad_vcf_url,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    DEFAULT_BATCH_ROWS,
    SOURCE_CONTIG_COL,
    PanelBatchContext,
    PanelWriteSummary,
    write_checked_panel,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import (
    GnomadGroup,
    panel_af_col,
    panel_an_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL
from mecfs_bio.util.subproc.run_command import execute_command_with_retries

logger = structlog.get_logger()

DEFAULT_MAX_ATTEMPTS = 6
GNOMAD_PART_SUB_FOLDER = "per_chromosome"
_TSV_NULL = "."
_INFO_ID = re.compile(r"^##INFO=<ID=([^,>]+)")
_CONTIG_ID = re.compile(r"^##contig=<ID=([^,>]+)")
_ASSEMBLY = re.compile(r"assembly=([^,>]+)")


@frozen(slots=True)
class GnomadChromosomeAlleleFrequencyTask(Task):
    meta: ReferenceFileMeta
    release: GnomadRelease
    chrom: int
    fasta_task: Task
    batch_rows: int = DEFAULT_BATCH_ROWS
    max_attempts: int = DEFAULT_MAX_ATTEMPTS

    @property
    def deps(self) -> list[Task]:
        return [self.fasta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        url = gnomad_vcf_url(self.release, self.chrom)
        fasta = load_fasta(fetch, self.fasta_task)
        assert_gnomad_header(
            read_vcf_header(url, self.max_attempts), self.release, self.chrom
        )
        tsv_path = scratch_dir / "records.tsv"
        execute_command_with_retries(
            gnomad_query_command(url, tsv_path, self.release),
            max_attempts=self.max_attempts,
        )
        out_path = scratch_dir / "allele_frequencies.parquet"
        summary = write_checked_panel(
            gnomad_tsv_batches(tsv_path, self.release, self.batch_rows),
            out_path,
            PanelBatchContext(
                fasta=fasta,
                contig_codes={gnomad_contig_name(self.release, self.chrom): self.chrom},
            ),
            byte_stream_split_columns=[
                panel_an_col(group) for group in self.release.groups
            ],
        )
        assert_gnomad_summary(summary, self.chrom)
        logger.info(
            "gnomAD chromosome allele frequencies written",
            release=self.release.name,
            chromosome=self.chrom,
            rows_written=summary.rows_written,
            fasta_ambiguous_dropped=summary.counts.fasta_ambiguous,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, release: GnomadRelease, chrom: int, fasta_task: Task
    ) -> "GnomadChromosomeAlleleFrequencyTask":
        fasta_meta = fasta_task.meta
        assert isinstance(fasta_meta, FASTAMeta), (
            f"fasta_task must carry FASTAMeta, got {type(fasta_meta).__name__}"
        )
        assert fasta_meta.build == release.build, (
            f"{release.name} is build {release.build}; the FASTA is build {fasta_meta.build}"
        )
        assert chrom in release.chromosomes, (
            f"chromosome {chrom} is not in {release.name}: {release.chromosomes}"
        )
        contig = gwaslab_code_to_contig_name(chrom)
        return cls(
            meta=ReferenceFileMeta(
                group=GNOMAD_GROUP,
                sub_group=release.name,
                sub_folder=PurePath(GNOMAD_PART_SUB_FOLDER),
                id=AssetId(f"gnomad_{release.name}_chr{contig}_allele_frequencies"),
                filename=f"chr{contig}_allele_frequencies",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            release=release,
            chrom=chrom,
            fasta_task=fasta_task,
        )


def read_vcf_header(url: str, max_attempts: int) -> str:
    return execute_command_with_retries(
        ["bcftools", "view", "-h", shlex.quote(url)], max_attempts=max_attempts
    )


def assert_gnomad_header(header: str, release: GnomadRelease, chrom: int) -> None:
    lines = header.splitlines()
    info_ids = {match.group(1) for line in lines if (match := _INFO_ID.match(line))}
    required = [
        column
        for group in release.groups
        for column in (panel_af_col(group), panel_an_col(group))
    ]
    missing = [column for column in required if column not in info_ids]
    assert not missing, f"{release.name} header lacks INFO fields {missing}"
    contig_lines = [line for line in lines if _CONTIG_ID.match(line)]
    assemblies = {
        match.group(1) for line in contig_lines if (match := _ASSEMBLY.search(line))
    }
    assert assemblies == {release.header_assembly}, (
        f"{release.name} header contig assemblies {sorted(assemblies)} != "
        f"{release.header_assembly!r}"
    )
    contigs = {
        match.group(1) for line in contig_lines if (match := _CONTIG_ID.match(line))
    }
    expected_contig = gnomad_contig_name(release, chrom)
    assert expected_contig in contigs, (
        f"{release.name} header does not declare contig {expected_contig}"
    )


def gnomad_query_command(
    url: str, tsv_path: Path, release: GnomadRelease
) -> list[str]:
    """One bcftools query process: PASS records polymorphic in a main group, as a TSV."""
    polymorphic = " || ".join(f"{panel_af_col(group)}>0" for group in release.main_groups)
    fields = "".join(
        f"\\t%INFO/{panel_af_col(group)}\\t%INFO/{panel_an_col(group)}"
        for group in release.groups
    )
    return [
        "bcftools",
        "query",
        "-i",
        f"'FILTER=\"PASS\" && ({polymorphic})'",
        "-f",
        f"'%CHROM\\t%POS\\t%REF\\t%ALT{fields}\\n'",
        "-o",
        shlex.quote(str(tsv_path)),
        shlex.quote(url),
    ]


def gnomad_tsv_batches(
    tsv_path: Path, release: GnomadRelease, batch_rows: int
) -> Iterator[pl.DataFrame]:
    """Stream the headerless bcftools TSV in batches of batch_rows rows."""
    assert tsv_path.stat().st_size > 0, (
        f"{release.name}: no PASS records polymorphic in a main group in {tsv_path}"
    )
    schema: dict[str, pl.DataType] = {
        SOURCE_CONTIG_COL: pl.String(),
        GWASLAB_POS_COL: pl.Int32(),
        PANEL_REF_COL: pl.String(),
        PANEL_ALT_COL: pl.String(),
    }
    for group in release.groups:
        schema[panel_af_col(group)] = pl.Float32()
        schema[panel_an_col(group)] = pl.Int32()
    records = pl.scan_csv(
        tsv_path,
        separator="\t",
        has_header=False,
        schema=schema,
        null_values=[_TSV_NULL],
        quote_char=None,
    )
    # collect_batches is marked unstable in polars; the Task tests pin its batching.
    for batch in records.collect_batches(chunk_size=batch_rows):
        _assert_af_null_exactly_where_an_is_zero(batch, release.groups)
        yield batch


def _assert_af_null_exactly_where_an_is_zero(
    frame: pl.DataFrame, groups: tuple[GnomadGroup, ...]
) -> None:
    for group in groups:
        af = pl.col(panel_af_col(group))
        an = pl.col(panel_an_col(group))
        bad = frame.filter(an.is_null() | (af.is_null() != (an == 0)))
        assert bad.height == 0, (
            f"{group}: AF must be null exactly where AN is 0, and AN never null; "
            f"e.g. {bad.head(3).rows()}"
        )


def assert_gnomad_summary(summary: PanelWriteSummary, chrom: int) -> None:
    assert summary.counts.ref_mismatch == 0, (
        f"chromosome {chrom}: {summary.counts.ref_mismatch} gnomAD REF alleles differ "
        "from a pure A/C/G/T FASTA span"
    )
    assert summary.rows_written > 0, f"chromosome {chrom}: the table is empty"
    assert summary.chromosomes == frozenset({chrom}), (
        f"chromosome {chrom}: rows on {sorted(summary.chromosomes)}"
    )
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_chromosome_allele_frequency_task.py -v`
Expected: all PASS.

- [ ] **Step 6: Run green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_task4.log
git add mecfs_bio/build_system/task/genome_reference_harmonization/gnomad test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_chromosome_allele_frequency_task.py
git commit -m "Add gnomAD release descriptor and per-chromosome allele-frequency Task

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: gnomAD panel concatenation Task and its task generator

**Files:**
- Create: `mecfs_bio/build_system/task/genome_reference_harmonization/gnomad/gnomad_allele_frequency_panel_task.py`
- Create: `mecfs_bio/build_system/task_generator/gnomad_allele_frequency_panel_task_generator.py`
- Test: `test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_allele_frequency_panel_task.py`

**Interfaces:**
- Consumes: `GnomadRelease`, `GNOMAD_GROUP`, `GnomadChromosomeAlleleFrequencyTask.create` (Task 4); `fetch_file_path` (Task 3); `parquet_encoding`, `open_parquet_writer` (Task 2); `PanelAlleleFrequencyColumns.prefixed`, `panel_an_col` (Task 1).
- Produces:
  - `GnomadAlleleFrequencyPanelTask.create(asset_id: str, release: GnomadRelease, part_tasks: Sequence[GnomadChromosomeAlleleFrequencyTask])`; fields `meta`, `release`, `part_tasks: tuple[GnomadChromosomeAlleleFrequencyTask, ...]`. The parts are injected, not built by create(); construction asserts they are exactly the release's chromosomes, in release order, from the same release.
  - `generate_gnomad_allele_frequency_panel_tasks(asset_id: str, release: GnomadRelease, fasta_task: Task) -> GnomadAlleleFrequencyPanelTasks`; `GnomadAlleleFrequencyPanelTasks` (frozen) has `part_tasks`, `panel_task`, `terminal_tasks()`. This is the repo's usual home for building a group of related Tasks: `build_system/task_generator` (generic, may use Task classes but no concrete asset instances), as opposed to `asset_generator` (may reference concrete assets).

- [ ] **Step 1: Write the failing tests**

`test_gnomad_allele_frequency_panel_task.py`:

```python
"""Task-level tests of GnomadAlleleFrequencyPanelTask on pre-written per-chromosome parts."""

from pathlib import Path, PurePath

import polars as pl
import pyarrow.parquet as pq
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_allele_frequency_panel_task import (
    GnomadAlleleFrequencyPanelTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    panel_af_col,
    panel_an_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

_RELEASE = GnomadRelease(
    name="test_release",
    build="19",
    vcf_url_template="unused/sites.{chrom}.vcf.bgz",
    contig_prefix="",
    chromosomes=(1, 23),
    main_groups=("afr", "nfe"),
    extra_groups=(),
    header_assembly="gnomAD_GRCh37",
)
_FASTA = FakeTask(
    FASTAMeta(
        group="genome_sequence",
        sub_group="synthetic",
        sub_folder=PurePath("processed"),
        id=AssetId("fasta"),
        build="19",
    )
)


def _part(chrom: int, positions: list[int], extra_column: bool = False) -> pl.DataFrame:
    n = len(positions)
    frame = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [chrom] * n,
            GWASLAB_POS_COL: positions,
            "REF": ["A"] * n,
            "ALT": ["G"] * n,
            panel_af_col("afr"): [0.1] * n,
            panel_an_col("afr"): [100] * n,
            panel_af_col("nfe"): [0.2] * n,
            panel_an_col("nfe"): [200] * n,
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            "REF": pl.String,
            "ALT": pl.String,
            panel_af_col("afr"): pl.Float32,
            panel_an_col("afr"): pl.Int32,
            panel_af_col("nfe"): pl.Float32,
            panel_an_col("nfe"): pl.Int32,
        },
    )
    if extra_column:
        frame = frame.with_columns(pl.lit(0).alias("UNEXPECTED"))
    return frame


def _part_task(chrom: int) -> GnomadChromosomeAlleleFrequencyTask:
    return GnomadChromosomeAlleleFrequencyTask.create(
        release=_RELEASE, chrom=chrom, fasta_task=_FASTA
    )


def _panel_task() -> GnomadAlleleFrequencyPanelTask:
    return GnomadAlleleFrequencyPanelTask.create(
        asset_id="gnomad_panel",
        release=_RELEASE,
        part_tasks=[_part_task(chrom) for chrom in _RELEASE.chromosomes],
    )


def _run(tmp_path: Path, parts: dict[int, pl.DataFrame]) -> Path:
    task = _panel_task()
    assets: dict[str, Asset] = {}
    for part_task in task.part_tasks:
        path = tmp_path / f"{part_task.asset_id}.parquet"
        parts[part_task.chrom].write_parquet(path)
        assets[part_task.asset_id] = FileAsset(path)

    def fetch(asset_id: AssetId) -> Asset:
        return assets[asset_id]

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return result.path


def test_parts_are_concatenated_in_chromosome_order(tmp_path: Path) -> None:
    path = _run(tmp_path, {23: _part(23, [5]), 1: _part(1, [3, 7])})
    panel = pl.read_parquet(path)
    assert panel.select(GWASLAB_CHROM_COL, GWASLAB_POS_COL).rows() == [
        (1, 3),
        (1, 7),
        (23, 5),
    ]
    columns = _panel_task().meta.allele_frequency_columns
    assert columns is not None
    assert columns.column_for("nfe") in panel.columns
    row_group = pq.ParquetFile(path).metadata.row_group(0)
    encodings = {
        row_group.column(i).path_in_schema: set(row_group.column(i).encodings)
        for i in range(row_group.num_columns)
    }
    assert "BYTE_STREAM_SPLIT" in encodings[panel_an_col("nfe")]
    assert "RLE_DICTIONARY" in encodings[panel_af_col("nfe")]


def test_parts_with_different_schemas_fail(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        _run(tmp_path, {1: _part(1, [3]), 23: _part(23, [5], extra_column=True)})


@pytest.mark.parametrize("chromosomes", [[1], [23, 1], [1, 23, 23]])
def test_parts_not_matching_the_release_chromosomes_are_rejected(
    chromosomes: list[int],
) -> None:
    # A missing, reordered or repeated part would give an incomplete or unsorted panel.
    with pytest.raises(AssertionError):
        GnomadAlleleFrequencyPanelTask.create(
            asset_id="gnomad_panel",
            release=_RELEASE,
            part_tasks=[_part_task(chrom) for chrom in chromosomes],
        )
```

- [ ] **Step 2: Run them to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_allele_frequency_panel_task.py -v`
Expected: collection error, ModuleNotFoundError.

- [ ] **Step 3: Implement**

`gnomad_allele_frequency_panel_task.py`:

```python
"""
A gnomAD release's allele-frequency panel: the per-chromosome tables concatenated in gwaslab
chromosome order into one parquet sorted by (CHR, POS), so the harmonizer's per-chromosome
filter prunes row groups.

The per-chromosome Tasks are injected by the caller (normally
generate_gnomad_allele_frequency_panel_tasks) and are ordinary dependencies, kept in the asset
store under their own sub_folder. They cannot be deleted to save space, because the build system materializes
every transitive dependency of a target; a path_remap rule can move them to another disk.
Do not wrap this Task in DiscardDepsWrapper: the multi-hour build would become
all-or-nothing again and the FASTA would be rebuilt in a temporary store.
"""

from pathlib import Path, PurePath
from typing import Sequence

import pyarrow
import pyarrow.parquet
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    open_parquet_writer,
    parquet_encoding,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GNOMAD_GROUP,
    GnomadRelease,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    fetch_file_path,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import panel_an_col
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL

GNOMAD_PANEL_FILENAME = "gnomad_allele_frequencies"


@frozen(slots=True)
class GnomadAlleleFrequencyPanelTask(Task):
    meta: HarmonizableReferenceTableMeta
    release: GnomadRelease
    part_tasks: tuple[GnomadChromosomeAlleleFrequencyTask, ...]

    def __attrs_post_init__(self) -> None:
        # Exactly the release's chromosomes, in release order: concatenation relies on it
        # for a complete panel sorted by (CHR, POS).
        part_chromosomes = tuple(part.chrom for part in self.part_tasks)
        assert part_chromosomes == self.release.chromosomes, (
            f"{self.release.name}: parts cover chromosomes {part_chromosomes}, "
            f"expected {self.release.chromosomes} in that order"
        )
        foreign = [
            part.asset_id for part in self.part_tasks if part.release != self.release
        ]
        assert not foreign, f"parts from a release other than {self.release.name}: {foreign}"

    @property
    def deps(self) -> list[Task]:
        return list(self.part_tasks)

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        part_paths = [fetch_file_path(fetch, part) for part in self.part_tasks]
        schema = identical_part_schema(part_paths)
        out_path = scratch_dir / (GNOMAD_PANEL_FILENAME + ".parquet")
        concatenate_parts(
            part_paths,
            out_path,
            schema,
            byte_stream_split_columns=[
                panel_an_col(group) for group in self.release.groups
            ],
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        release: GnomadRelease,
        part_tasks: Sequence[GnomadChromosomeAlleleFrequencyTask],
    ) -> "GnomadAlleleFrequencyPanelTask":
        return cls(
            meta=HarmonizableReferenceTableMeta(
                group=GNOMAD_GROUP,
                sub_group=release.name,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename=GNOMAD_PANEL_FILENAME,
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=release.build,
                    ref_allele_col=PANEL_REF_COL,
                    pos_col=GWASLAB_POS_COL,
                ),
                allele_frequency_columns=PanelAlleleFrequencyColumns.prefixed(
                    release.groups
                ),
            ),
            release=release,
            part_tasks=tuple(part_tasks),
        )


def identical_part_schema(part_paths: list[Path]) -> pyarrow.Schema:
    """The parts' shared arrow schema (metadata ignored); fails if any part differs."""
    schemas = [pyarrow.parquet.read_schema(path).remove_metadata() for path in part_paths]
    assert schemas, "no per-chromosome parts"
    differing = [
        str(path) for path, schema in zip(part_paths, schemas) if not schema.equals(schemas[0])
    ]
    assert not differing, f"parts whose schema differs from {part_paths[0]}: {differing}"
    return schemas[0]


def concatenate_parts(
    part_paths: list[Path],
    out_path: Path,
    schema: pyarrow.Schema,
    byte_stream_split_columns: list[str],
) -> None:
    """Copy every row group of every part, in order, into one parquet."""
    encoding = parquet_encoding(
        schema.names,
        compression="zstd",
        compression_level=None,
        byte_stream_split_columns=byte_stream_split_columns,
    )
    with open_parquet_writer(out_path, schema, encoding) as writer:
        for path in part_paths:
            part = pyarrow.parquet.ParquetFile(path)
            for index in range(part.num_row_groups):
                writer.write_table(part.read_row_group(index).replace_schema_metadata(None))
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_allele_frequency_panel_task.py -v`
Expected: PASS. (Parts written by polars carry int32/float32 and large_string; the writer's schema is the parts' schema, so the tables match it.)

- [ ] **Step 5: Add the task generator**

`mecfs_bio/build_system/task_generator/gnomad_allele_frequency_panel_task_generator.py`. It is wiring only, so it has no test of its own (repo convention); Task 6's asset module constructs it at import time.

```python
"""
Builds a gnomAD release's allele-frequency panel together with its per-chromosome parts.
"""

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_allele_frequency_panel_task import (
    GnomadAlleleFrequencyPanelTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    GnomadChromosomeAlleleFrequencyTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)


@frozen(slots=True)
class GnomadAlleleFrequencyPanelTasks:
    part_tasks: tuple[GnomadChromosomeAlleleFrequencyTask, ...]
    panel_task: GnomadAlleleFrequencyPanelTask

    def terminal_tasks(self) -> list[Task]:
        return [self.panel_task]


def generate_gnomad_allele_frequency_panel_tasks(
    asset_id: str, release: GnomadRelease, fasta_task: Task
) -> GnomadAlleleFrequencyPanelTasks:
    part_tasks = tuple(
        GnomadChromosomeAlleleFrequencyTask.create(
            release=release, chrom=chrom, fasta_task=fasta_task
        )
        for chrom in release.chromosomes
    )
    return GnomadAlleleFrequencyPanelTasks(
        part_tasks=part_tasks,
        panel_task=GnomadAlleleFrequencyPanelTask.create(
            asset_id=asset_id, release=release, part_tasks=part_tasks
        ),
    )
```

- [ ] **Step 6: Run green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_task5.log
git add mecfs_bio/build_system/task/genome_reference_harmonization/gnomad mecfs_bio/build_system/task_generator/gnomad_allele_frequency_panel_task_generator.py test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_gnomad_allele_frequency_panel_task.py
git commit -m "Add gnomAD allele-frequency panel concatenation Task

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Assets and a live header probe

**Files:**
- Create: `mecfs_bio/assets/reference_data/gnomad/__init__.py` (empty), `gnomad_releases.py`, `gnomad_allele_frequency_panels.py`
- Create: `mecfs_bio/assets/reference_data/pan_ukbb/__init__.py` (empty), `pan_ukbb_variant_manifest.py`, `pan_ukbb_allele_frequencies.py`
- Create: `experiments/claude/gnomad_af_reference/check_release_headers.py`

**Interfaces:**
- Consumes: `GnomadRelease`, `generate_gnomad_allele_frequency_panel_tasks`, `read_vcf_header`, `assert_gnomad_header`, `gnomad_vcf_url` (Tasks 4-5); `PanUkbbAlleleFrequencyPanelTask.create` (Task 3); `DownloadFileTask`; `UCSC_HG19_INDEXED_FASTA`, `UCSC_HG38_INDEXED_FASTA`.
- Produces: `GNOMAD_V2_1_1_GENOMES`, `GNOMAD_V4_1_GENOMES`, `GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES`, `GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES`, `PAN_UKBB_VARIANT_MANIFEST`, `PAN_UKBB_HG19_ALLELE_FREQUENCIES`.

There is no unit test for asset wiring (repo convention: import-time construction, which green's import of all asset modules exercises, is enough). The live probe in Step 4 is the check.

- [ ] **Step 1: gnomAD releases**

`mecfs_bio/assets/reference_data/gnomad/gnomad_releases.py`:

```python
"""
gnomAD genome sites-VCF releases used for allele-frequency panels.

Files are in the public bucket gnomad-public-us-east-1 (anonymous access, no egress charge
to us). Groups and header facts were measured on chr21; see
experiments/claude/design_specs/2026-10-02-gnomad-allele-frequency-panel-design.md.
"""

from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GnomadRelease,
)

_AUTOSOMES = tuple(range(1, 23))
_X = 23
_Y = 24

GNOMAD_V2_1_1_GENOMES = GnomadRelease(
    name="v2_1_1_genomes",
    build="19",
    vcf_url_template=(
        "https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/vcf/genomes/"
        "gnomad.genomes.r2.1.1.sites.{chrom}.vcf.bgz"
    ),
    contig_prefix="",
    chromosomes=(*_AUTOSOMES, _X),
    main_groups=("afr", "amr", "asj", "eas", "fin", "nfe"),
    extra_groups=("oth", "nfe_nwe", "nfe_seu", "nfe_onf", "nfe_est"),
    header_assembly="gnomAD_GRCh37",
)

GNOMAD_V4_1_GENOMES = GnomadRelease(
    name="v4_1_genomes",
    build="38",
    vcf_url_template=(
        "https://gnomad-public-us-east-1.s3.amazonaws.com/release/4.1/vcf/genomes/"
        "gnomad.genomes.v4.1.sites.chr{chrom}.vcf.bgz"
    ),
    contig_prefix="chr",
    chromosomes=(*_AUTOSOMES, _X, _Y),
    main_groups=("afr", "amr", "asj", "eas", "fin", "nfe", "sas"),
    extra_groups=("ami", "mid", "remaining"),
    header_assembly="gnomAD_GRCh38",
)
```

- [ ] **Step 2: gnomAD panels**

`gnomad_allele_frequency_panels.py`:

```python
"""
gnomAD genome allele-frequency panels for genome-reference harmonization.

The hg19 panel is built from v2.1.1 (about 14 h streaming, about 3.7 GiB). The hg38 panel
from v4.1 is defined but built only on request (about 14 h, about 10 GiB). Neither is wrapped
in DiscardDepsWrapper; the per-chromosome parts stay in the asset store under
reference_data/gnomad/<release>/per_chromosome/, which a path_remap rule may move.
"""

from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_releases import (
    GNOMAD_V2_1_1_GENOMES,
    GNOMAD_V4_1_GENOMES,
)
from mecfs_bio.build_system.task_generator.gnomad_allele_frequency_panel_task_generator import (
    generate_gnomad_allele_frequency_panel_tasks,
)

GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES = (
    generate_gnomad_allele_frequency_panel_tasks(
        asset_id="gnomad_v2_1_1_genomes_hg19_allele_frequencies",
        release=GNOMAD_V2_1_1_GENOMES,
        fasta_task=UCSC_HG19_INDEXED_FASTA,
    ).panel_task
)

GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES = (
    generate_gnomad_allele_frequency_panel_tasks(
        asset_id="gnomad_v4_1_genomes_hg38_allele_frequencies",
        release=GNOMAD_V4_1_GENOMES,
        fasta_task=UCSC_HG38_INDEXED_FASTA,
    ).panel_task
)
```

- [ ] **Step 3: Pan-UKBB assets**

`mecfs_bio/assets/reference_data/pan_ukbb/pan_ukbb_variant_manifest.py`:

```python
"""
The Pan-UK Biobank variant manifest (full_variant_qc_metrics.txt.bgz, 2.7 GB, GRCh37),
pinned by md5 (equal to its S3 ETag; last modified 2020-08-28).

It stays in the asset store as a dependency of the Pan-UKBB panel; its own sub_folder
(raw/) lets a path_remap rule move it to another disk.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    PAN_UKBB_MANIFEST_READ_SPEC,
)

PAN_UKBB_VARIANT_MANIFEST = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="pan_ukbb_variant_manifest",
        group="pan_ukbb",
        sub_group="variant_manifest",
        sub_folder=PurePath("raw"),
        extension=".txt.bgz",
        read_spec=PAN_UKBB_MANIFEST_READ_SPEC,
    ),
    url=(
        "https://pan-ukb-us-east-1.s3.amazonaws.com/sumstats_release/"
        "full_variant_qc_metrics.txt.bgz"
    ),
    md5_hash="e70ebc8289f762dd8d5086f54e766654",
)
```

`pan_ukbb_allele_frequencies.py`:

```python
"""
Pan-UK Biobank allele-frequency panel (GRCh37, imputed UK Biobank, six genetic-ancestry
groups) for genome-reference harmonization.

295 manifest rows (chr21 47, chr22 86, X 162) have ref and alt swapped relative to the hg19
FASTA, inherited from UK Biobank's imputed BGEN allele order; the panel drops them and the
build fails if the count changes. Measured by experiments/claude/pan_ukbb_manifest/.
"""

from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_variant_manifest import (
    PAN_UKBB_VARIANT_MANIFEST,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.pan_ukbb.pan_ukbb_allele_frequency_panel_task import (
    PanUkbbAlleleFrequencyPanelTask,
)

PAN_UKBB_HG19_ALLELE_FREQUENCIES = PanUkbbAlleleFrequencyPanelTask.create(
    asset_id="pan_ukbb_hg19_allele_frequencies",
    manifest_task=PAN_UKBB_VARIANT_MANIFEST,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    groups=("ukb_afr", "ukb_amr", "ukb_csa", "ukb_eas", "ukb_eur", "ukb_mid"),
    chromosomes=(*range(1, 23), 23),
    expected_ref_mismatches=295,
)
```

- [ ] **Step 4: Live header probe**

`experiments/claude/gnomad_af_reference/check_release_headers.py`:

```python
"""
Check every per-chromosome URL of both gnomAD releases against its release description.

Runs the gnomAD Task's own header check (INFO fields for every group, contig assembly,
contig declared) on each chromosome's live header, so a wrong URL template or contig name
fails here in minutes rather than hours into a build.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.check_release_headers \
        2>&1 | tee experiments/claude/gnomad_af_reference/check_release_headers.log
"""

from mecfs_bio.assets.reference_data.gnomad.gnomad_releases import (
    GNOMAD_V2_1_1_GENOMES,
    GNOMAD_V4_1_GENOMES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_chromosome_allele_frequency_task import (
    assert_gnomad_header,
    read_vcf_header,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    gnomad_vcf_url,
)


def main() -> None:
    for release in [GNOMAD_V2_1_1_GENOMES, GNOMAD_V4_1_GENOMES]:
        for chrom in release.chromosomes:
            url = gnomad_vcf_url(release, chrom)
            assert_gnomad_header(read_vcf_header(url, max_attempts=3), release, chrom)
            print(f"ok {release.name} chromosome {chrom} {url}")


if __name__ == "__main__":
    main()
```

Run: `pixi r python -m experiments.claude.gnomad_af_reference.check_release_headers 2>&1 | tee experiments/claude/gnomad_af_reference/check_release_headers.log`
Expected: 47 "ok" lines (23 for v2.1.1, 24 for v4.1). If one fails, fix the release description in gnomad_releases.py, not the check.

- [ ] **Step 5: Run green and commit**

```bash
pixi r invoke green 2>&1 | tee /tmp/claude-1000/green_task6.log
git add mecfs_bio/assets/reference_data/gnomad mecfs_bio/assets/reference_data/pan_ukbb experiments/claude/gnomad_af_reference/check_release_headers.py
git commit -m "Add gnomAD and Pan-UKBB allele-frequency panel assets

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Validation on real data

These steps build real assets in the shared asset store, run from the main checkout with its existing `default_runner_config.yaml`. **Before Step 1, ask the user** whether another build job is using the store: two processes writing the same info store at once can corrupt it.

**Files:**
- Create: `experiments/claude/gnomad_af_reference/byte_identity_1000g.py`, `build_panel.py`, `decode_me_panel_comparison.py`
- Modify: the spec (a Validation results section), and the memory note project_pan_ukbb_variant_manifest_af_source.md

- [ ] **Step 1: Byte identity of an existing 1000 Genomes harmonization**

`byte_identity_1000g.py`:

```python
"""
Validation step 1: the panel-ancestry refactor leaves 1000 Genomes harmonizations unchanged.

Reads DecodeME's current harmonized table (built by the pre-refactor code), force-rebuilds it
with the current code, and compares the two as frames and as bytes. The build cache does not
track code, hence the forced rebuild.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.byte_identity_1000g \
        2>&1 | tee experiments/claude/gnomad_af_reference/byte_identity_1000g.log
"""

import hashlib
from pathlib import Path

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task


def _built_path(task: Task, rebuild: bool) -> Path:
    assets = DEFAULT_RUNNER.run(
        [task], must_rebuild_transitive=[task] if rebuild else ()
    )
    asset = assets[task.asset_id]
    assert isinstance(asset, FileAsset)
    return asset.path


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    task = DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task
    before_path = _built_path(task, rebuild=False)
    before_frame = pl.read_parquet(before_path)
    before_hash = _sha256(before_path)
    after_path = _built_path(task, rebuild=True)
    print(f"rows before {before_frame.height}")
    print("frames equal:", before_frame.equals(pl.read_parquet(after_path)))
    print("bytes equal:", before_hash == _sha256(after_path))


if __name__ == "__main__":
    main()
```

Run it as in the docstring. Expected: "frames equal: True". "bytes equal" should also be True; if only bytes differ, report it to the user (parquet writer metadata) rather than treating it as a failure.

- [ ] **Step 2: Build the panels**

`build_panel.py`:

```python
"""
Validation step 2: build one allele-frequency panel and report its size and layout.

Usage (one panel per run):
    pixi r python -m experiments.claude.gnomad_af_reference.build_panel pan_ukbb \
        2>&1 | tee experiments/claude/gnomad_af_reference/build_panel_pan_ukbb.log
    pixi r python -m experiments.claude.gnomad_af_reference.build_panel gnomad_v2 \
        2>&1 | tee experiments/claude/gnomad_af_reference/build_panel_gnomad_v2.log
The gnomAD v2.1.1 build streams about 451 GiB and takes about 14 h; it resumes per
chromosome if interrupted. The Task logs (rows written, drops) are in the tee'd log.
"""

import sys

import polars as pl
import pyarrow.parquet as pq

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_allele_frequencies import (
    PAN_UKBB_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

PANELS: dict[str, Task] = {
    "pan_ukbb": PAN_UKBB_HG19_ALLELE_FREQUENCIES,
    "gnomad_v2": GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
}


def main(name: str) -> None:
    task = PANELS[name]
    asset = DEFAULT_RUNNER.run([task])[task.asset_id]
    assert isinstance(asset, FileAsset)
    metadata = pq.ParquetFile(asset.path).metadata
    print(f"{name}: {asset.path}")
    print(f"size {asset.path.stat().st_size / 2**30:.2f} GiB, rows {metadata.num_rows}, "
          f"row groups {metadata.num_row_groups}")
    print(
        pl.scan_parquet(asset.path)
        .group_by(GWASLAB_CHROM_COL)
        .len()
        .sort(GWASLAB_CHROM_COL)
        .collect()
    )


if __name__ == "__main__":
    assert len(sys.argv) == 2 and sys.argv[1] in PANELS, f"usage: build_panel.py {list(PANELS)}"
    main(sys.argv[1])
```

Run pan_ukbb first (minutes). Expected: about 28,987,239 rows, 23 chromosomes, the Task's log line showing ref_mismatches_dropped=295. Then start gnomad_v2 in the background (run_in_background, timeout up to 7,200,000 ms; re-launch if the timeout ends it, since parts already built are reused). Expected at the end: about 210 M rows, about 3.7 GiB, fasta_ambiguous counts in the per-chromosome log lines.

- [ ] **Step 3: DecodeME four-way comparison**

`decode_me_panel_comparison.py`:

```python
"""
Validation step 3: harmonize DecodeME against four panel/ancestry choices and compare.

Choices: 1000 Genomes "eur", gnomAD v2.1.1 "nfe" and "nfe_nwe", Pan-UKBB "ukb_eur". For each,
the production trust decision and per-chromosome resolution run on DecodeME's
pre-harmonization table. Reported: the trust decision, drop-reason counts (kept, NOT_IN_PANEL,
AF_MISMATCH, ...), per-variant agreement of the outcome (kept orientation or drop reason)
between choices, and palindromes kept or flipped on a chosen-ancestry AF of exactly 0.

Usage:
    pixi r python -m experiments.claude.gnomad_af_reference.decode_me_panel_comparison \
        2>&1 | tee experiments/claude/gnomad_af_reference/decode_me_panel_comparison.log
"""

import itertools

import polars as pl
from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.pan_ukbb.pan_ukbb_allele_frequencies import (
    PAN_UKBB_HG19_ALLELE_FREQUENCIES,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    PanelTable,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    resolve_chromosome_rows,
    resolve_panel_af_col,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    decide_trust,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.allele_frequency_panel_constants import PanelAncestry
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
OUTCOME_COL = "outcome"
KEPT = "kept"
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
# The DecodeME rsID-assignment chain runs with palindromes dropped when unresolved.
OPTIONS = GenomeReferenceHarmonizationOptions()


@frozen(slots=True)
class PanelChoice:
    label: str
    panel_task: Task
    ancestry: PanelAncestry


CHOICES = [
    PanelChoice("1000g_eur", THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES, "eur"),
    PanelChoice("gnomad_nfe", GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES, "nfe"),
    PanelChoice("gnomad_nfe_nwe", GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES, "nfe_nwe"),
    PanelChoice("pan_ukbb_eur", PAN_UKBB_HG19_ALLELE_FREQUENCIES, "ukb_eur"),
]


def _file(asset: Asset) -> FileAsset:
    assert isinstance(asset, FileAsset)
    return asset


def resolve_all(
    choice: PanelChoice, sumstats: pl.LazyFrame, fasta: IndexedFasta, panel: PanelTable
) -> pl.DataFrame:
    """Every input row, keyed by SNPID, with its outcome: oriented EA/NEA or a drop reason."""
    rules = resolve_column_rules(columns=sumstats.collect_schema().names(), extra=())
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
    evidence = count_trust_evidence_genome_wide(
        sumstats, chromosomes, fasta, panel, OPTIONS
    )
    trusted = decide_trust(evidence, OPTIONS)
    print(f"[{choice.label}] trusted={trusted} suspicious={evidence.suspicious}")
    parts = []
    for chrom in chromosomes:
        context = ChromosomeContext(
            chrom=chrom, fasta=fasta, trusted=trusted, rules=rules, options=OPTIONS
        )
        resolved = resolve_chromosome_rows(sumstats, context, panel)
        parts.append(
            resolved.select(
                GWASLAB_SNPID_COL,
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                EA,
                NEA,
                pl.coalesce(
                    pl.col(DROP_REASON_COL),
                    pl.lit(KEPT + ":") + pl.col(EA) + pl.lit("/") + pl.col(NEA),
                ).alias(OUTCOME_COL),
            )
        )
    return pl.concat(parts)


def report_drop_reasons(label: str, outcomes: pl.DataFrame) -> None:
    reason = (
        pl.when(pl.col(OUTCOME_COL).str.starts_with(KEPT))
        .then(pl.lit(KEPT))
        .otherwise(pl.col(OUTCOME_COL))
    )
    print(f"\n== [{label}] outcome counts")
    with pl.Config(tbl_rows=30):
        print(outcomes.group_by(reason.alias("reason")).len().sort("len", descending=True))


def report_agreement(outcomes: dict[str, pl.DataFrame]) -> None:
    print("\n== pairwise outcome disagreement (rows whose outcome differs)")
    for left, right in itertools.combinations(outcomes, 2):
        joined = outcomes[left].select(GWASLAB_SNPID_COL, OUTCOME_COL).join(
            outcomes[right].select(GWASLAB_SNPID_COL, OUTCOME_COL),
            on=GWASLAB_SNPID_COL,
            suffix="_right",
        )
        differing = joined.filter(pl.col(OUTCOME_COL) != pl.col(OUTCOME_COL + "_right"))
        print(f"{left} vs {right}: {differing.height} of {joined.height}")
        with pl.Config(tbl_rows=10, fmt_str_lengths=40):
            print(
                differing.group_by(OUTCOME_COL, OUTCOME_COL + "_right")
                .len()
                .sort("len", descending=True)
                .head(10)
            )


def report_zero_af_palindromes(
    label: str, outcomes: pl.DataFrame, panel: PanelTable
) -> None:
    palindromic = (pl.col(EA).str.len_bytes() == 1) & (
        pl.col(EA).replace_strict(COMPLEMENT, default=None) == pl.col(NEA)
    )
    resolved = outcomes.filter(palindromic & pl.col(OUTCOME_COL).str.starts_with(KEPT))
    panel_rows = pl.scan_parquet(panel.path).select(
        GWASLAB_CHROM_COL,
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
        PANEL_REF_COL,
        PANEL_ALT_COL,
        pl.col(panel.af_col).alias("panel_af"),
    )
    matched = (
        resolved.lazy()
        .with_columns(pl.col(GWASLAB_CHROM_COL).cast(pl.Int32), pl.col(GWASLAB_POS_COL).cast(pl.Int64))
        .join(panel_rows, on=[GWASLAB_CHROM_COL, GWASLAB_POS_COL])
        .filter(
            ((pl.col(PANEL_REF_COL) == pl.col(NEA)) & (pl.col(PANEL_ALT_COL) == pl.col(EA)))
            | ((pl.col(PANEL_REF_COL) == pl.col(EA)) & (pl.col(PANEL_ALT_COL) == pl.col(NEA)))
        )
        .collect(engine="streaming")
    )
    print(
        f"\n== [{label}] resolved palindromes: {resolved.height}; "
        f"decided by a panel AF of exactly 0: {matched.filter(pl.col('panel_af') == 0).height}"
    )


def main() -> None:
    pre_task = DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.pre_harmonization_table_task
    panel_tasks = list({choice.panel_task.asset_id: choice.panel_task for choice in CHOICES}.values())
    assets = DEFAULT_RUNNER.run([pre_task, UCSC_HG19_INDEXED_FASTA, *panel_tasks])
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(assets[pre_task.asset_id], pre_task.meta, IdentityPipe())
    outcomes: dict[str, pl.DataFrame] = {}
    for choice in CHOICES:
        panel = PanelTable(
            path=_file(assets[choice.panel_task.asset_id]).path,
            af_col=resolve_panel_af_col(choice.panel_task, choice.ancestry),
        )
        outcomes[choice.label] = resolve_all(choice, sumstats, fasta, panel)
        report_drop_reasons(choice.label, outcomes[choice.label])
        report_zero_af_palindromes(choice.label, outcomes[choice.label], panel)
    report_agreement(outcomes)


if __name__ == "__main__":
    main()
```

Run it as in the docstring (after Step 2's gnomAD build finishes). Expected: four outcome tables, six pairwise disagreement lines, four zero-AF palindrome counts. Check that GWASLAB_SNPID_COL exists in DecodeME's pre-harmonization table first (`pl.read_parquet_schema` on the asset); if it does not, key the agreement join on (CHR, POS, input EA, input NEA) captured before resolution instead.

- [ ] **Step 4: Record results**

Add a "## Validation results (2026-10-..)" section to the spec, after "Validation on real data", with: byte-identity outcome; Pan-UKBB size and row count; gnomAD v2.1.1 size, row count and total FASTA-ambiguous drops; the four outcome tables (kept / NOT_IN_PANEL / AF_MISMATCH and other reasons) side by side; the pairwise disagreement counts; zero-AF palindrome counts; and a short paragraph on what they imply for the deferred AF = 0 palindrome question and for choosing a default panel. Do not edit docs/. Update the memory note project_pan_ukbb_variant_manifest_af_source.md with the built panel's size.

- [ ] **Step 5: Commit**

```bash
git add experiments/claude/gnomad_af_reference/*.py experiments/claude/design_specs/2026-10-02-gnomad-allele-frequency-panel-design.md
git commit -m "Validate gnomAD and Pan-UKBB panels on DecodeME

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
