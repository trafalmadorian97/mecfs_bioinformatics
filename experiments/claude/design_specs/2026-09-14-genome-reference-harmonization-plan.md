# Genome-Reference Harmonization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace gwaslab harmonization with an in-repo, streaming, polars
Task. The Task orients summary statistics against a genome FASTA and a
reference panel, and never keeps an indel whose orientation it cannot
establish.

**Architecture:**

- **Reference inputs.** Two pinned reference assets feed the Task: an indexed
  uncompressed FASTA, and a sites-only allele-frequency parquet built from the
  panel VCF.
- **Library.** A pure-function library under
  mecfs_bio/build_system/task/genome_reference_harmonization/ does the work:
  1. memory-mapped FASTA gather;
  2. allele classes;
  3. the trust decision;
  4. palindrome and ambiguous-indel rules;
  5. the column flip registry.
- **Task.** GenomeReferenceHarmonizationTask streams the input one chromosome
  at a time in two passes:
  - pass 1: the trust decision;
  - pass 2: resolution, writing one parquet part per chromosome and
    concatenating them with a streaming sink.

**Tech Stack:**

- Python 3.13, polars 1.44, narwhals 2.26, numpy memmap, pysam (faidx),
  bcftools (pixi), attrs, structlog.
- pytest with typeguard, run via pixi.

**Spec:** experiments/claude/design_specs/2026-09-14-genome-reference-harmonization-design.md.
Read it before starting; this plan argues from it.

## Global Constraints

**Running things**

- Run every command through pixi: pixi r python ..., pixi r invoke green.
  Capture green output to a logfile and grep the pytest summary; never pipe it
  to tail.
- Work on a new branch genome-reference-harmonization created from main. The
  current branch susie-test-debug has an unrelated uncommitted test edit;
  leave it alone.
- Commit only when the user has approved executing this plan.

**Coding conventions (from project memory)**

- **Docstrings:** no backticks around inline code, no RST.
- **Types:**
  - use Path objects, not strings, except at subprocess and JSON boundaries;
  - Literal aliases for enumerable string values, each value also exposed as a Final constant typed with its alias. Code uses the constants, never the raw strings, and an import-time assertion checks alias and constants agree. A StrEnum is not used, because polars turns Python Enum members into its own Enum dtype inside expressions.
  - column names are constants too, including test-only and experiment labels;
  - no Any, no Callable[...]; use Protocols;
  - return frozen attrs objects, never bare tuples.
- **Structure:**
  - module-level helper free functions, not helper methods;
  - pass same-typed arguments by keyword;
  - column names come from constants (GWASLAB_*_COL in
    mecfs_bio/constants/gwaslab_constants.py, or module constants defined
    here), never repeated literals.
- **Fail fast:** assert preconditions with clear messages. Cross-field option
  invariants go in __attrs_post_init__.
- **Dataframes:** polars for new code. Subprocesses go through execute_command
  (mecfs_bio/util/subproc/run_command.py).
- **Task.create():**
  - derives trait and project from an input task's meta;
  - production Tasks never use SimpleFileMeta;
  - reference data uses ReferenceFileMeta or ReferenceDataDirectoryMeta.
- **Tests:**
  - Task-level unit tests: FakeTask dependencies, synthetic inputs in
    tmp_path, a fetch mapping asset id to Asset, task.execute(...,
    wf=make_wf());
  - no monkeypatching or mocks; inject dependencies;
  - never assert on error-message or log text: pytest.raises(AssertionError)
    without match=;
  - no cargo-cult tests; define shared constants instead of retyping literals;
  - no structural tests for asset wiring modules (import-time construction is
    enough);
  - no skipif guards on library presence.
- **Terminology** in code and prose:
  - "gwaslab harmonization" (the old step);
  - "genome-reference harmonization" (this Task);
  - "task-reference harmonization" (HarmonizeGWASWithReferenceViaAlleles,
    untouched).

**Values copied from the spec**

- Trust bar: exactly 100%. That is, consistent_snvs >= min_checkable_snvs
  (default 10,000) and zero inconsistent SNVs, and consistent_indels >=
  min_checkable_indels (default 1,000) and zero inconsistent indels.
- Palindrome thresholds: sumstats MAF 0.4 and panel MAF 0.4, epsilon 1e-6
  (gwaslab defaults).
- Stringent indel options: indel_max_af_distance 0.1 and indel_min_af_margin
  0.2 until Task 9 (experiment V3) sets final values.
- excluded_chromosomes default: (25,) (MT).
- Pinned md5s:
  - hg19.fa.gz 806c02398f5ac5da8ffd6da2d1d5d1a9 (UCSC md5sum.txt);
  - hg38.fa.gz 1c9dcaddfa41027f17cd8f7a82c7293b (UCSC md5sum.txt);
  - 1kg EUR hg19 VCF 2c78cb84cb1f90b576510decc45e5b9b (gwaslab catalogue).

**Clarifications to the spec (decided while planning; mention in the PR)**

1. **Invalid alleles are dropped, not asserted.** Rows whose uppercased EA or
   NEA is not ^[ACGT]+$, or has EA == NEA, are dropped with reason
   invalid_allele. gwaslab basic_check does not remove such rows, so asserting
   would fail on real inputs. Null CHR, POS, EA or NEA are still asserted.
2. **Parsimony is not asserted.** gwaslab normalization leaves complex indels
   (e.g. AC to GTT) untrimmed. They are classified like other
   different-length variants.
3. **Panel coverage is not asserted.** A chromosome with no panel rows (e.g.
   Y in the hg19 panel) simply has no records: untrusted palindromes and
   ambiguous indels there are dropped. FASTA coverage of every non-excluded
   chromosome is asserted.
4. **indel_min_af_margin must be strictly positive**, so the keep and flip
   readings can never both be chosen.
5. **Output dtypes:** POS is cast to Int64 and EA/NEA to String. STATUS is
   dropped from the output.

## File Structure

**Create: library** (mecfs_bio/build_system/task/genome_reference_harmonization/)

| File | Responsibility |
|---|---|
| __init__.py | empty |
| fasta.py | FaiEntry, IndexedFasta, contig_to_gwaslab_code, reference_matches (memory-mapped gather) |
| indexed_fasta_task.py | IndexedFastaTask: gz FASTA to directory with genome.fa and .fai |
| reference_panel_task.py | ReferencePanelAlleleFrequencyTask: VCF to sites-only AF parquet |
| flip.py | FlipRule, BoundPair, ExtraColumnRule, COLUMN_FLIP_RULES, resolve_column_rules, flip_statistics |
| options.py | GenomeReferenceHarmonizationOptions |
| allele_classes.py | AlleleClass, reverse_complement_expr, prepare_alleles, valid_alleles_expr, classify_alleles |
| trust.py | TrustCounts, count_trust_evidence, decide_trust |
| ambiguous_indels.py | decide_ambiguous_indels (stringent rules) |
| palindromes.py | decide_palindrome_strands |
| outcomes.py | AlleleAction, PalindromeDecision, DropReason Literal aliases and their ACTION_*, PALINDROME_*, DROP_* constants |
| resolve_chromosome.py | PanelLoader, ChromosomeContext, resolve_chromosome |
| genome_reference_harmonization_task.py | GenomeReferenceHarmonizationTask, scan_sumstats_as_polars, chromosomes_to_harmonize, count_trust_evidence_genome_wide, ParquetPanelLoader |

**Create: pipe, assets, experiments**

| File | Responsibility |
|---|---|
| mecfs_bio/build_system/task/pipes/drop_indels_pipe.py | DropIndelsPipe |
| mecfs_bio/assets/reference_data/genome_sequence/__init__.py | empty |
| .../genome_sequence/ucsc_hg19_fasta.py | UCSC_HG19_FASTA_GZ, UCSC_HG19_INDEXED_FASTA |
| .../genome_sequence/ucsc_hg38_fasta.py | UCSC_HG38_FASTA_GZ, UCSC_HG38_INDEXED_FASTA |
| mecfs_bio/assets/reference_data/thousand_genomes/eur_hg19_vcf.py | THOUSAND_GENOMES_EUR_HG19_VCF |
| .../thousand_genomes/eur_panel_allele_frequencies.py | THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES, THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES |
| experiments/claude/genome_reference_harmonization/ | V1-V5 scripts and logs |

**Create: tests** (test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/)

| File | Responsibility |
|---|---|
| __init__.py | empty |
| genome_reference_fixtures.py | synthetic genome, panel writer, Task runner |
| test_indexed_fasta_task.py | IndexedFastaTask |
| test_reference_panel_task.py | ReferencePanelAlleleFrequencyTask |
| test_genome_reference_harmonization_task.py | the Task |
| test_mecfs_bio/unit/build_system/task/pipes/test_drop_indels_pipe.py | DropIndelsPipe (create the pipes test dir with __init__.py if absent) |

**Modify**

- mecfs_bio/asset_generator/annovar_37_basic_rsid_assignment.py
- The Liu IBD harmonized asset, its dump asset (deleted) and the two consumers
  of the dump asset.
- test_mecfs_bio/system/test_harmonize_drop_ambiguous.py and
  .github/workflows/system_test_harmonize.yml
- mecfs_bio/build_system/task/gwaslab/gwaslab_create_sumstats_task.py (remove
  dead harmonization code, Task 12)

---

### Task 1: Memory-mapped FASTA reader and IndexedFastaTask

**Files:**
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/__init__.py (empty)
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/fasta.py
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/indexed_fasta_task.py
- Create: test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/__init__.py (empty)
- Create: test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_indexed_fasta_task.py

**Interfaces:**
- Produces:
  - FASTA_FILENAME = "genome.fa"
  - DEFAULT_MAX_GATHER_BYTES: int
  - FaiEntry(length, offset, line_bases, line_bytes)
  - contig_to_gwaslab_code(name: str) -> int | None
  - IndexedFasta(fasta_path: Path, entries: Mapping[int, FaiEntry]), with classmethod IndexedFasta.open(directory: Path)
  - reference_matches(fasta: IndexedFasta, chrom: int, positions: np.ndarray, alleles: pl.Series, max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES) -> np.ndarray (bool)
  - IndexedFastaTask.create(fasta_gz_task: Task, asset_id: str) -> IndexedFastaTask, producing a DirectoryAsset

- [ ] **Step 1: Write the failing test**

```python
"""Task-level test of IndexedFastaTask, exercising the memory-mapped reference gather."""

import gzip
from pathlib import Path, PurePath

import numpy as np
import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.indexed_fasta_task import (
    IndexedFastaTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf

_LINE_WIDTH = 7
# 1-based: position 1 is A; 31-40 is a run of ten T.
_CHR1 = "ACGTTGGGGGCATGCAATTCGATCGATCGATTTTTTTTTTACACACACACGTCAGTCAGT"
_CHRX = "ggggaaaacccc"  # lowercase (soft-masked) bases must still match
_GZ_ID = "fasta_gz"


def _write_gzipped_fasta(path: Path) -> None:
    lines: list[str] = []
    for name, sequence in {"chr1": _CHR1, "chrX": _CHRX, "chrUn_gl000220": "ACGT"}.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[start : start + _LINE_WIDTH]
            for start in range(0, len(sequence), _LINE_WIDTH)
        )
    with gzip.open(path, "wt") as out:
        out.write("\n".join(lines) + "\n")


def _build_indexed_fasta(tmp_path: Path) -> IndexedFasta:
    gz_path = tmp_path / "genome.fa.gz"
    _write_gzipped_fasta(gz_path)
    source = FakeTask(
        ReferenceFileMeta(
            group="genome_sequence",
            sub_group="synthetic",
            sub_folder=PurePath("raw"),
            extension=".fa.gz",
            id=AssetId(_GZ_ID),
        )
    )
    task = IndexedFastaTask.create(fasta_gz_task=source, asset_id="indexed")

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == _GZ_ID
        return FileAsset(gz_path)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    return IndexedFasta.open(result.path)


def test_indexed_fasta_maps_main_contigs_to_gwaslab_codes(tmp_path: Path) -> None:
    fasta = _build_indexed_fasta(tmp_path)
    assert set(fasta.entries) == {1, 23}


@pytest.mark.parametrize("max_gather_bytes", [1, 1_000_000])
def test_reference_matches_across_line_wraps(tmp_path: Path, max_gather_bytes: int) -> None:
    fasta = _build_indexed_fasta(tmp_path)
    positions = np.array([1, 5, 5, 31, 31, 58, 60], dtype=np.int64)
    alleles = pl.Series(["A", "TG", "TC", "T" * 10, "T" * 11, "AGT", "TA"])
    matched = reference_matches(
        fasta,
        chrom=1,
        positions=positions,
        alleles=alleles,
        max_gather_bytes=max_gather_bytes,
    )
    # TA at 60 runs past the contig end, so it cannot match.
    assert matched.tolist() == [True, True, False, True, False, True, False]


def test_reference_matches_is_case_insensitive(tmp_path: Path) -> None:
    fasta = _build_indexed_fasta(tmp_path)
    matched = reference_matches(
        fasta,
        chrom=23,
        positions=np.array([1, 5], dtype=np.int64),
        alleles=pl.Series(["GGGG", "aaaac"]),
    )
    assert matched.tolist() == [True, True]
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_indexed_fasta_task.py -q`
Expected: FAIL with ModuleNotFoundError for genome_reference_harmonization.fasta.

- [ ] **Step 3: Implement fasta.py**

```python
"""
Memory-mapped access to an uncompressed, faidx-indexed genome FASTA.

The uncompressed file is memory-mapped, and (chromosome, position) is translated to a
byte offset with the .fai index, so nothing is loaded up front. The operating system
page cache serves the touched pages, and those pages are reclaimable rather than
anonymous memory. Alleles are compared as prefixes of the reference starting at the
variant position, ignoring case (soft-masked bases are lowercase).
"""

import re
from collections.abc import Mapping
from pathlib import Path

import numpy as np
import polars as pl
from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME

FASTA_FILENAME = "genome.fa"
FAI_SUFFIX = ".fai"
# Upper bound on the byte size of one gather's offset matrix, so that very long alleles
# are processed in row slices instead of materializing a huge matrix.
DEFAULT_MAX_GATHER_BYTES = 64 * 1024 * 1024
_OFFSET_BYTES = 8
_ASCII_UPPERCASE_MASK = 0xDF
_MAIN_CONTIG = re.compile(r"^(?:chr)?([0-9]+|X|Y|M|MT)$")


@frozen
class FaiEntry:
    length: int # Number of actual nucleotide bases
    offset: int # Location in fasta file where first base starts
    line_bases: int # Each normal line contains this many bases
    line_bytes: int # Each normal line contains this many bytes


def contig_to_gwaslab_code(name: str) -> int | None:
    """gwaslab numeric chromosome code for a main contig name (chr1, 1, chrX, chrM, MT), else None."""
    match = _MAIN_CONTIG.match(name)
    if match is None:
        return None
    token = match.group(1)
    if token == "M":
        token = "MT"
    if token in GWASLAB_CHROM_CODE_FOR_NAME:
        return GWASLAB_CHROM_CODE_FOR_NAME[token]
    return int(token)


@frozen
class IndexedFasta:
    """An uncompressed FASTA and its .fai entries, keyed by gwaslab chromosome code."""

    fasta_path: Path
    entries: Mapping[int, FaiEntry]

    @classmethod
    def open(cls, directory: Path) -> "IndexedFasta":
        fasta_path = directory / FASTA_FILENAME
        fai_path = directory / (FASTA_FILENAME + FAI_SUFFIX)
        assert fasta_path.is_file(), f"missing {fasta_path}"
        assert fai_path.is_file(), f"missing {fai_path}"
        entries: dict[int, FaiEntry] = {}
        for line in fai_path.read_text().splitlines():
            name, length, offset, line_bases, line_bytes = line.split("\t")[:5]
            code = contig_to_gwaslab_code(name)
            if code is None:
                continue
            assert code not in entries, f"two FASTA contigs map to chromosome {code}"
            entries[code] = FaiEntry(
                length=int(length),
                offset=int(offset),
                line_bases=int(line_bases),
                line_bytes=int(line_bytes),
            )
        return cls(fasta_path=fasta_path, entries=entries)


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
    order = np.argsort(lengths, kind="stable")
    boundaries = np.flatnonzero(np.diff(lengths[order])) + 1
    for group in np.split(order, boundaries):
        length = int(lengths[group[0]])
        rows_per_slice = max(1, max_gather_bytes // (length * _OFFSET_BYTES))
        for start in range(0, len(group), rows_per_slice):
            rows = group[start : start + rows_per_slice]
            result[rows] = _match_equal_length(
                genome=genome,
                entry=entry,
                positions=positions[rows],
                alleles=alleles.gather(rows),
                length=length,
            )
    return result


def _match_equal_length(
    genome: np.ndarray,
    entry: FaiEntry,
    positions: np.ndarray,
    alleles: pl.Series,
    length: int,
) -> np.ndarray:
    start = positions.astype(np.int64) - 1
    in_bounds = (start >= 0) & (start + length <= entry.length)
    clipped = np.clip(start, 0, max(entry.length - length, 0))
    base_index = clipped[:, None] + np.arange(length, dtype=np.int64)[None, :]
    file_offsets = (
        entry.offset
        + base_index
        + (base_index // entry.line_bases) * (entry.line_bytes - entry.line_bases)
    )
    reference = genome[file_offsets] & _ASCII_UPPERCASE_MASK
    observed = (
        np.frombuffer("".join(alleles.to_list()).encode("ascii"), dtype=np.uint8).reshape(
            len(positions), length
        )
        & _ASCII_UPPERCASE_MASK
    )
    return in_bounds & (reference == observed).all(axis=1)
```

- [ ] **Step 4: Implement indexed_fasta_task.py**

```python
"""
Decompress a gzipped genome FASTA and index it with faidx.

Memory-mapped reference lookups need the uncompressed file plus its .fai, so the
output is a directory holding both. Wrap instances in DiscardDepsWrapper so that the
gzipped download is not also kept in the asset store.
"""

import gzip
import shutil
from pathlib import Path, PurePath

import pysam
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    FASTA_FILENAME,
)
from mecfs_bio.build_system.wf.base_wf import WF

_COPY_CHUNK_BYTES = 16 * 1024 * 1024


@frozen
class IndexedFastaTask(Task):
    meta: ReferenceDataDirectoryMeta
    fasta_gz_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.fasta_gz_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.fasta_gz_task.asset_id)
        assert isinstance(source, FileAsset), (
            f"expected {self.fasta_gz_task.asset_id} to be a FileAsset, got {type(source).__name__}"
        )
        out_dir = scratch_dir / "indexed_fasta"
        out_dir.mkdir()
        fasta_path = out_dir / FASTA_FILENAME
        with gzip.open(source.path, "rb") as compressed, open(fasta_path, "wb") as plain:
            shutil.copyfileobj(compressed, plain, length=_COPY_CHUNK_BYTES)
        pysam.faidx(str(fasta_path))
        return DirectoryAsset(out_dir)

    @classmethod
    def create(cls, fasta_gz_task: Task, asset_id: str) -> "IndexedFastaTask":
        source_meta = fasta_gz_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ReferenceDataDirectoryMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
            ),
            fasta_gz_task=fasta_gz_task,
        )
```

- [ ] **Step 5: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_indexed_fasta_task.py -q`
Expected: 4 passed.

- [ ] **Step 6: Commit**

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization test_mecfs_bio/unit/build_system/task/genome_reference_harmonization
git commit -m "Add memory-mapped FASTA reader and IndexedFastaTask"
```

### Task 2: ReferencePanelAlleleFrequencyTask

**Files:**
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/reference_panel_task.py
- Create: test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_reference_panel_task.py

**Interfaces:**
- Consumes: contig_to_gwaslab_code from Task 1.
- Produces:
  - PANEL_REF_COL = "REF", PANEL_ALT_COL = "ALT", PANEL_AF_COL = "AF"
  - PANEL_COLUMNS = [CHR, POS, REF, ALT, AF]
  - Panel parquet schema: CHR Int32, POS Int32, REF String, ALT String, AF Float32, sorted by CHR then POS, one row per (CHR, POS, REF, ALT)
  - ReferencePanelAlleleFrequencyTask.create(vcf_task: Task, asset_id: str)

- [ ] **Step 1: Write the failing test**

```python
"""Task-level test of ReferencePanelAlleleFrequencyTask on a tiny VCF."""

from pathlib import Path, PurePath

import polars as pl
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_COLUMNS,
    PANEL_REF_COL,
    ReferencePanelAlleleFrequencyTask,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL

_VCF_ID = "panel_vcf"
_VCF = (
    "##fileformat=VCFv4.2\n"
    '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n'
    "##contig=<ID=1>\n##contig=<ID=X>\n##contig=<ID=GL000191.1>\n"
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
    "1\t200\t.\tTG\tT\t.\t.\tAF=0.86\tGT\t0|1\n"
    "1\t100\t.\tA\tG\t.\t.\tAF=0.2\tGT\t0|0\n"
    "1\t200\t.\tT\tTG\t.\t.\tAF=0\tGT\t0|0\n"  # mirrored record at the same site: kept
    "1\t300\t.\tC\tT\t.\t.\tAF=0.1\tGT\t0|0\n"
    "1\t300\t.\tC\tT\t.\t.\tAF=0.1\tGT\t0|0\n"  # exact duplicate: collapsed
    "1\t400\t.\tG\tA\t.\t.\tAF=0.3\tGT\t0|0\n"
    "1\t400\t.\tG\tA\t.\t.\tAF=0.4\tGT\t0|0\n"  # conflicting duplicate: removed
    "1\t500\t.\tA\t.\t.\t.\tAF=.\tGT\t0|0\n"  # monomorphic, no ALT: removed
    "X\t50\t.\tC\tG\t.\t.\tAF=0.5\tGT\t0|1\n"
    "GL000191.1\t10\t.\tA\tC\t.\t.\tAF=0.5\tGT\t0|1\n"  # non-main contig: removed
)


def _run(tmp_path: Path) -> pl.DataFrame:
    vcf_path = tmp_path / "panel.vcf"
    vcf_path.write_text(_VCF)
    source = FakeTask(
        ReferenceFileMeta(
            group="thousand_genomes",
            sub_group="synthetic",
            sub_folder=PurePath("raw"),
            extension=".vcf",
            id=AssetId(_VCF_ID),
        )
    )
    task = ReferencePanelAlleleFrequencyTask.create(vcf_task=source, asset_id="panel_af")

    def fetch(asset_id: AssetId) -> Asset:
        assert asset_id == _VCF_ID
        return FileAsset(vcf_path)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def test_panel_keeps_main_contig_sites_with_unambiguous_frequencies(tmp_path: Path) -> None:
    panel = _run(tmp_path)
    assert panel.columns == PANEL_COLUMNS
    assert panel.rows() == [
        (1, 100, "A", "G", pytest.approx(0.2)),
        (1, 200, "T", "TG", pytest.approx(0.0)),
        (1, 200, "TG", "T", pytest.approx(0.86)),
        (1, 300, "C", "T", pytest.approx(0.1)),
        (23, 50, "C", "G", pytest.approx(0.5)),
    ]
    assert panel.schema[GWASLAB_CHROM_COL] == pl.Int32
    assert panel.schema[GWASLAB_POS_COL] == pl.Int32
    assert panel.schema[PANEL_AF_COL] == pl.Float32
    assert panel.schema[PANEL_REF_COL] == pl.String
    assert panel.schema[PANEL_ALT_COL] == pl.String
```

The two mirrored rows at position 200 compare by REF within POS, so the sort key must be (CHR, POS, REF, ALT) for a deterministic order.

- [ ] **Step 2: Run the test to verify it fails**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_reference_panel_task.py -q`
Expected: FAIL with ModuleNotFoundError for reference_panel_task.

- [ ] **Step 3: Implement reference_panel_task.py**

```python
"""
Convert a reference panel VCF into a sites-only allele-frequency parquet.

Genome-reference harmonization looks up panel allele frequencies for palindromic SNVs
and for indels whose alleles both match the genome. This Task strips genotypes with
bcftools and keeps CHR (gwaslab numeric coding), POS, REF, ALT and AF.

- Only main contigs are kept (1-22, X, Y, MT).
- Records with no ALT or no AF are dropped.
- The panel is expected to be split and normalized already. A multi-allelic record's
  comma-separated AF fails the Float32 parse, which stops the build.
- Exact duplicate (CHR, POS, REF, ALT) keys with one AF collapse to one row. Keys with
  conflicting AF are removed and counted in the log, so a lookup never sees two answers.

Rows are sorted by CHR, POS so a per-chromosome scan reads only matching row groups.
The work is done one contig at a time to bound memory. Wrap instances in
DiscardDepsWrapper so the multi-gigabyte VCF is not also kept in the asset store.
bcftools is resolved from the pixi environment.
"""

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
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    contig_to_gwaslab_code,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

PANEL_REF_COL = "REF"
PANEL_ALT_COL = "ALT"
PANEL_AF_COL = "AF"
PANEL_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    PANEL_REF_COL,
    PANEL_ALT_COL,
    PANEL_AF_COL,
]
_CONTIG_COL = "contig"
_N_AF_COL = "n_distinct_af"
_SITE_KEYS = [GWASLAB_POS_COL, PANEL_REF_COL, PANEL_ALT_COL]


def _write_sites_tsv(vcf_path: Path, tsv_path: Path) -> None:
    execute_command(
        [
            "bcftools",
            "view",
            "-G",
            "-Ou",
            str(vcf_path),
            "|",
            "bcftools",
            "query",
            "-f",
            r"'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n'",
            "-o",
            str(tsv_path),
        ]
    )


def _scan_sites(tsv_path: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        tsv_path,
        separator="\t",
        has_header=False,
        null_values=["."],
        schema={
            _CONTIG_COL: pl.String,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            PANEL_AF_COL: pl.Float32,
        },
    )


def _contigs_by_code(sites: pl.LazyFrame) -> dict[int, list[str]]:
    names = sites.select(pl.col(_CONTIG_COL).unique()).collect(engine="streaming")[_CONTIG_COL]
    by_code: dict[int, list[str]] = {}
    skipped: list[str] = []
    for name in names.to_list():
        code = contig_to_gwaslab_code(name)
        if code is None:
            skipped.append(name)
        else:
            by_code.setdefault(code, []).append(name)
    if skipped:
        logger.info("skipping non-main panel contigs", contigs=sorted(skipped))
    return by_code


def _one_chromosome(sites: pl.LazyFrame, code: int, names: list[str]) -> pl.DataFrame:
    rows = (
        sites.filter(pl.col(_CONTIG_COL).is_in(names))
        .drop_nulls([PANEL_ALT_COL, PANEL_AF_COL])
        .collect(engine="streaming")
    )
    grouped = rows.group_by(_SITE_KEYS).agg(
        pl.col(PANEL_AF_COL).first(),
        pl.col(PANEL_AF_COL).n_unique().alias(_N_AF_COL),
    )
    n_conflicting = grouped.filter(pl.col(_N_AF_COL) > 1).height
    if n_conflicting:
        logger.warning(
            "removing panel sites with conflicting duplicate allele frequencies",
            chromosome=code,
            n_sites=n_conflicting,
        )
    return (
        grouped.filter(pl.col(_N_AF_COL) == 1)
        .with_columns(pl.lit(code, dtype=pl.Int32).alias(GWASLAB_CHROM_COL))
        .select(PANEL_COLUMNS)
        .sort(_SITE_KEYS)
    )


@frozen
class ReferencePanelAlleleFrequencyTask(Task):
    meta: ReferenceFileMeta
    vcf_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.vcf_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        vcf = fetch(self.vcf_task.asset_id)
        assert isinstance(vcf, FileAsset), (
            f"expected {self.vcf_task.asset_id} to be a FileAsset, got {type(vcf).__name__}"
        )
        tsv_path = scratch_dir / "panel_sites.tsv"
        _write_sites_tsv(vcf_path=vcf.path, tsv_path=tsv_path)
        sites = _scan_sites(tsv_path)
        parts_dir = scratch_dir / "parts"
        parts_dir.mkdir()
        part_paths: list[Path] = []
        for code, names in sorted(_contigs_by_code(sites).items()):
            part_path = parts_dir / f"chr{code}.parquet"
            _one_chromosome(sites, code, names).write_parquet(part_path)
            part_paths.append(part_path)
        assert part_paths, f"no main-contig records in {vcf.path}"
        out_path = scratch_dir / "panel_allele_frequencies.parquet"
        pl.concat([pl.scan_parquet(path) for path in part_paths]).sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(cls, vcf_task: Task, asset_id: str) -> "ReferencePanelAlleleFrequencyTask":
        source_meta = vcf_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=ReferenceFileMeta(
                group="reference_panel_allele_frequencies",
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="panel_allele_frequencies",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            vcf_task=vcf_task,
        )
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_reference_panel_task.py -q`
Expected: 1 passed. If bcftools rejects the VCF's AF=0 or AF=. lines, fix the test VCF, not the Task.

- [ ] **Step 5: Commit**

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization/reference_panel_task.py test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_reference_panel_task.py
git commit -m "Add ReferencePanelAlleleFrequencyTask"
```

### Task 3: Pinned reference assets (hg19 and hg38 FASTA; 1kg EUR panels)

These are wiring-only asset modules. Per project convention there are no unit
tests for them; import-time construction plus one real build is the check.

**Files:**
- Create: mecfs_bio/assets/reference_data/genome_sequence/__init__.py (empty)
- Create: mecfs_bio/assets/reference_data/genome_sequence/ucsc_hg19_fasta.py
- Create: mecfs_bio/assets/reference_data/genome_sequence/ucsc_hg38_fasta.py
- Create: mecfs_bio/assets/reference_data/thousand_genomes/eur_hg19_vcf.py
- Create: mecfs_bio/assets/reference_data/thousand_genomes/eur_panel_allele_frequencies.py
- Create: experiments/claude/genome_reference_harmonization/build_reference_assets.py

**Interfaces:**
- Consumes: IndexedFastaTask (Task 1), ReferencePanelAlleleFrequencyTask (Task 2), the existing THOUSAND_GENOMES_EUR_HG38_30X_VCF.
- Produces:
  - UCSC_HG19_INDEXED_FASTA and UCSC_HG38_INDEXED_FASTA (DiscardDepsWrapper; DirectoryAsset)
  - THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES and THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES (DiscardDepsWrapper; FileAsset)

- [ ] **Step 1: Write ucsc_hg19_fasta.py**

```python
"""
UCSC hg19 (GRCh37) genome sequence, pinned by the md5 UCSC publishes in md5sum.txt.

The indexed form is what genome-reference harmonization reads. It is wrapped in
DiscardDepsWrapper, so only the uncompressed FASTA and its .fai are stored, not the
gzipped download as well. Note that UCSC hg19 chrM is not the rCRS mitochondrial
sequence, which is why harmonization excludes MT by default.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.genome_reference_harmonization.indexed_fasta_task import (
    IndexedFastaTask,
)

UCSC_HG19_FASTA_GZ = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="ucsc_hg19_fasta_gz",
        group="genome_sequence",
        sub_group="ucsc_hg19",
        sub_folder=PurePath("raw"),
        extension=".fa.gz",
    ),
    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz",
    md5_hash="806c02398f5ac5da8ffd6da2d1d5d1a9",
)

UCSC_HG19_INDEXED_FASTA = DiscardDepsWrapper(
    IndexedFastaTask.create(
        fasta_gz_task=UCSC_HG19_FASTA_GZ, asset_id="ucsc_hg19_indexed_fasta"
    )
)
```

- [ ] **Step 2: Write ucsc_hg38_fasta.py**

```python
"""
UCSC hg38 (GRCh38) genome sequence, pinned by the md5 UCSC publishes in md5sum.txt.

Used to tune the stringent ambiguous-indel rules on build-38 DecodeME, whose source
orientation is fully reference-consistent. Wrapped in DiscardDepsWrapper so only the
uncompressed FASTA and its .fai are stored.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.genome_reference_harmonization.indexed_fasta_task import (
    IndexedFastaTask,
)

UCSC_HG38_FASTA_GZ = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="ucsc_hg38_fasta_gz",
        group="genome_sequence",
        sub_group="ucsc_hg38",
        sub_folder=PurePath("raw"),
        extension=".fa.gz",
    ),
    url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    md5_hash="1c9dcaddfa41027f17cd8f7a82c7293b",
)

UCSC_HG38_INDEXED_FASTA = DiscardDepsWrapper(
    IndexedFastaTask.create(
        fasta_gz_task=UCSC_HG38_FASTA_GZ, asset_id="ucsc_hg38_indexed_fasta"
    )
)
```

- [ ] **Step 3: Write eur_hg19_vcf.py**

```python
"""
The 1000 Genomes phase 3 EUR reference VCF (hg19), as distributed by gwaslab.

Multi-allelic variants are decomposed, variants are normalized, and INFO/AF carries
the EUR allele frequency. Records are indexed for chromosomes 1-22 and X only.

gwaslab republishes its panels in place (its hg38 panel silently gained chrX), so the
checksum is what pins the content. If the URL stops serving this checksum, rehost the
file as was done for eur_hg38_30x_vcf.py.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

THOUSAND_GENOMES_EUR_HG19_VCF = DownloadFileTask(
    meta=ReferenceFileMeta(
        id="thousand_genomes_eur_hg19_vcf",
        group="thousand_genomes",
        sub_group="eur_hg19",
        sub_folder=PurePath("raw"),
        extension=".vcf.gz",
    ),
    url=(
        "https://www.dropbox.com/scl/fi/sxnjd37t7e677wlluzeit/"
        "EUR.ALL.split_norm_af.1kgp3v5.hg19.vcf.gz"
        "?rlkey=9om2qu3tfjnq8nxj0tgwp8bd0&dl=1"
    ),
    md5_hash="2c78cb84cb1f90b576510decc45e5b9b",
)
```

- [ ] **Step 4: Write eur_panel_allele_frequencies.py**

```python
"""
Sites-only allele-frequency tables of the 1000 Genomes EUR panels, for
genome-reference harmonization.

Wrapped in DiscardDepsWrapper so the multi-gigabyte genotype VCFs are not kept in the
asset store, only the derived parquet.
"""

from mecfs_bio.assets.reference_data.thousand_genomes.eur_hg19_vcf import (
    THOUSAND_GENOMES_EUR_HG19_VCF,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_hg38_30x_vcf import (
    THOUSAND_GENOMES_EUR_HG38_30X_VCF,
)
from mecfs_bio.build_system.task.discard_deps_task_wrapper import DiscardDepsWrapper
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    ReferencePanelAlleleFrequencyTask,
)

THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES = DiscardDepsWrapper(
    ReferencePanelAlleleFrequencyTask.create(
        vcf_task=THOUSAND_GENOMES_EUR_HG19_VCF,
        asset_id="thousand_genomes_eur_hg19_panel_allele_frequencies",
    )
)

THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES = DiscardDepsWrapper(
    ReferencePanelAlleleFrequencyTask.create(
        vcf_task=THOUSAND_GENOMES_EUR_HG38_30X_VCF,
        asset_id="thousand_genomes_eur_hg38_panel_allele_frequencies",
    )
)
```

- [ ] **Step 5: Write the build script**

```python
"""Materialize the genome-reference harmonization reference assets and report their size.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.build_reference_assets \
    2>&1 | tee experiments/claude/genome_reference_harmonization/build_reference_assets.log
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

FASTAS = [UCSC_HG19_INDEXED_FASTA, UCSC_HG38_INDEXED_FASTA]
PANELS = [
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
]


def main() -> None:
    assets = DEFAULT_RUNNER.run([*FASTAS, *PANELS])
    for task in FASTAS:
        asset = assets[task.asset_id]
        assert isinstance(asset, DirectoryAsset)
        fasta = IndexedFasta.open(asset.path)
        print(task.asset_id, "chromosomes:", sorted(fasta.entries))
    for task in PANELS:
        asset = assets[task.asset_id]
        assert isinstance(asset, FileAsset)
        per_chrom = (
            pl.scan_parquet(asset.path)
            .group_by(GWASLAB_CHROM_COL)
            .len()
            .sort(GWASLAB_CHROM_COL)
            .collect()
        )
        print(task.asset_id, f"{asset.path.stat().st_size / 1e9:.2f} GB")
        print(per_chrom)


if __name__ == "__main__":
    main()
```

- [ ] **Step 6: Build the assets**

This downloads about 3.5 GB (hg19 panel), 1 GB (hg19 FASTA), 1 GB (hg38 FASTA)
and the hg38 panel if it is not already cached. It takes tens of minutes, so
run it in the background.

Run: `pixi r python -m experiments.claude.genome_reference_harmonization.build_reference_assets 2>&1 | tee experiments/claude/genome_reference_harmonization/build_reference_assets.log`

Expected:
- hg19 FASTA chromosomes [1..25]: the 22 autosomes plus X, Y, M.
- hg19 panel row counts for 1-22 and 23 only.
- hg38 panel including 23.
- Any "removing panel sites with conflicting duplicate allele frequencies"
  warnings: record their counts in the log.

- [ ] **Step 7: Commit**

```bash
git add mecfs_bio/assets/reference_data/genome_sequence mecfs_bio/assets/reference_data/thousand_genomes/eur_hg19_vcf.py mecfs_bio/assets/reference_data/thousand_genomes/eur_panel_allele_frequencies.py experiments/claude/genome_reference_harmonization
git commit -m "Add pinned FASTA and 1kg EUR panel allele-frequency reference assets"
```

### Task 4: Core GenomeReferenceHarmonizationTask

This covers allele classes, flips, trust counting and streaming execution.

**Temporary state after this task.** When a table is untrusted, palindromic
SNVs and ambiguous indels (class indel_both) are still kept in source
orientation, as in trusted mode:

- Task 5 adds untrusted palindrome strand rules.
- Task 6 adds the stringent ambiguous-indel rules.

Nothing wires the Task into an analysis until Task 11, so this is safe.

**Files:**
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/options.py
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/outcomes.py
- Modify: mecfs_bio/constants/gwaslab_constants.py (standard column-name constants)
- Modify: mecfs_bio/constants/regenie_constants.py (binary-trait column-name constants)
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/flip.py
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/allele_classes.py
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/trust.py
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/resolve_chromosome.py
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/genome_reference_harmonization_task.py
- Create: test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/genome_reference_fixtures.py
- Create: test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py

**Interfaces:**
- Consumes: IndexedFasta, reference_matches, FASTA_FILENAME, DEFAULT_MAX_GATHER_BYTES (Task 1); PANEL_* constants (Task 2).
- Produces:
  - GenomeReferenceHarmonizationOptions (all fields listed in options.py below)
  - FlipRule, InvertedBoundPair, ExtraColumnRule(column: str, rule: FlipRule), resolve_column_rules(columns, extra) -> Mapping[str, FlipRule], flip_statistics(frame, mask: pl.Expr, rules) -> pl.DataFrame, DROPPED_COLUMNS
  - ALLELE_CLASS_COL, IS_PALINDROMIC_SNV_COL, IS_PALINDROMIC_MNP_COL, AlleleClass, reverse_complement_expr(column), prepare_alleles(frame), valid_alleles_expr(), classify_alleles(frame, fasta, chrom, max_gather_bytes)
  - TrustCounts(consistent_snvs, inconsistent_snvs, consistent_indels, inconsistent_indels), count_trust_evidence(classified) -> TrustCounts, decide_trust(counts, options) -> bool
  - outcomes.py: AlleleAction, PalindromeDecision, DropReason and constants ACTION_KEEP, ACTION_SWAP, ACTION_COMPLEMENT, ACTION_COMPLEMENT_SWAP, PALINDROME_KEEP, PALINDROME_STRAND_FLIP, PALINDROME_UNRESOLVED, DROP_* (one per DropReason)
  - allele_classes.py constants CLASS_* (one per AlleleClass); flip.py constants FLIP_NEGATE, FLIP_COMPLEMENT, FLIP_INVERT, FLIP_INVARIANT
  - DROP_REASON_COL, ROW_INDEX_COL, PanelLoader, ChromosomeContext, resolve_chromosome(rows, context, load_panel) -> pl.DataFrame
  - GenomeReferenceHarmonizationTask.create(asset_id, sumstats_task, fasta_task, panel_task, options=..., pipe=...); scan_sumstats_as_polars, chromosomes_to_harmonize, count_trust_evidence_genome_wide, ParquetPanelLoader, resolve_chromosome_rows

- [ ] **Step 1: Write the test fixtures module**

```python
"""Synthetic genome, reference panel and Task runner for genome-reference harmonization tests."""

from collections.abc import Sequence
from pathlib import Path, PurePath

import polars as pl
import pysam
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    FASTA_FILENAME,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

# Line width 7 so that alleles cross FASTA line breaks. 1-based landmarks:
#   1-10  A C G T T G G G G G   (5: T followed by a G run -> T/TG is ambiguous)
#   11-20 C A T G C A A T T C   (13, 18, 19: T for A/T palindromes)
#   21-30 G A T C G A T C G A
#   31-40 T x 10                (a homopolymer for T/TT ambiguous indels)
#   41-50 A C A C A C A C A C   (41: AC for the AC/GT palindromic MNP)
#   51-60 G T C A G T C A G T
LINE_WIDTH = 7
CHR1_SEQUENCE = "ACGTTGGGGGCATGCAATTCGATCGATCGATTTTTTTTTTACACACACACGTCAGTCAGT"
CHR2_SEQUENCE = "GGGGAAAACCCCTTTT"
SE_VALUE = 0.01

TEST_OPTIONS = GenomeReferenceHarmonizationOptions(
    min_checkable_snvs=1, min_checkable_indels=1
)

_SUMSTATS_ID = "sumstats"
_FASTA_ID = "fasta"
_PANEL_ID = "panel"


@frozen
class Variant:
    pos: int
    ea: str
    nea: str
    eaf: float | None = 0.3
    beta: float = 0.1
    chrom: int = 1


@frozen
class PanelRecord:
    pos: int
    ref: str
    alt: str
    af: float
    chrom: int = 1


# A table of only these two is fully reference-consistent (trusted under TEST_OPTIONS);
# adding INCONSISTENT_SNV (EA is the reference base) makes it untrusted.
CONSISTENT_SNV = Variant(pos=1, ea="G", nea="A")
CONSISTENT_INDEL = Variant(pos=21, ea="GC", nea="G")
INCONSISTENT_SNV = Variant(pos=2, ea="C", nea="T", beta=0.2)


def sumstats_frame(variants: Sequence[Variant]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [v.chrom for v in variants],
            GWASLAB_POS_COL: [v.pos for v in variants],
            GWASLAB_EFFECT_ALLELE_COL: [v.ea for v in variants],
            GWASLAB_NON_EFFECT_ALLELE_COL: [v.nea for v in variants],
            GWASLAB_EFFECT_ALLELE_FREQ_COL: [v.eaf for v in variants],
            GWASLAB_BETA_COL: [v.beta for v in variants],
            GWASLAB_SE_COL: [SE_VALUE for _ in variants],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int64,
            GWASLAB_POS_COL: pl.Int64,
            GWASLAB_EFFECT_ALLELE_COL: pl.String,
            GWASLAB_NON_EFFECT_ALLELE_COL: pl.String,
            GWASLAB_EFFECT_ALLELE_FREQ_COL: pl.Float64,
            GWASLAB_BETA_COL: pl.Float64,
            GWASLAB_SE_COL: pl.Float64,
        },
    )


def _write_fasta(directory: Path) -> None:
    directory.mkdir()
    fasta_path = directory / FASTA_FILENAME
    lines: list[str] = []
    for name, sequence in {"chr1": CHR1_SEQUENCE, "chr2": CHR2_SEQUENCE}.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[start : start + LINE_WIDTH]
            for start in range(0, len(sequence), LINE_WIDTH)
        )
    fasta_path.write_text("\n".join(lines) + "\n")
    pysam.faidx(str(fasta_path))


def _write_panel(path: Path, records: Sequence[PanelRecord]) -> None:
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [r.chrom for r in records],
            GWASLAB_POS_COL: [r.pos for r in records],
            PANEL_REF_COL: [r.ref for r in records],
            PANEL_ALT_COL: [r.alt for r in records],
            PANEL_AF_COL: [r.af for r in records],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            PANEL_AF_COL: pl.Float32,
        },
    ).sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL).write_parquet(path)


def run_harmonization(
    work_dir: Path,
    sumstats: pl.DataFrame,
    panel: Sequence[PanelRecord] = (),
    options: GenomeReferenceHarmonizationOptions = TEST_OPTIONS,
    pipe: DataProcessingPipe = IdentityPipe(),
) -> pl.DataFrame:
    """Execute the Task on synthetic inputs in a fresh work_dir and return the output table."""
    work_dir.mkdir(parents=True)
    sumstats_path = work_dir / "sumstats.parquet"
    sumstats.write_parquet(sumstats_path)
    fasta_dir = work_dir / "fasta"
    _write_fasta(fasta_dir)
    panel_path = work_dir / "panel.parquet"
    _write_panel(panel_path, panel)
    parquet_spec = DataFrameReadSpec(DataFrameParquetFormat())
    task = GenomeReferenceHarmonizationTask.create(
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
            ReferenceDataDirectoryMeta(
                group="genome_sequence",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                id=AssetId(_FASTA_ID),
            )
        ),
        panel_task=FakeTask(
            ReferenceFileMeta(
                group="reference_panel_allele_frequencies",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                extension=".parquet",
                id=AssetId(_PANEL_ID),
                read_spec=parquet_spec,
            )
        ),
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


def positions(frame: pl.DataFrame, chrom: int = 1) -> list[int]:
    return frame.filter(pl.col(GWASLAB_CHROM_COL) == chrom)[GWASLAB_POS_COL].to_list()


def row_at(frame: pl.DataFrame, pos: int, chrom: int = 1) -> dict[str, object]:
    matching = frame.filter(
        (pl.col(GWASLAB_CHROM_COL) == chrom) & (pl.col(GWASLAB_POS_COL) == pos)
    )
    assert matching.height == 1, f"expected one row at {chrom}:{pos}, got {matching.height}"
    return matching.row(0, named=True)
```

- [ ] **Step 2: Write the failing Task tests**

```python
"""Task-level tests of GenomeReferenceHarmonizationTask on a synthetic genome and panel."""

from pathlib import Path

import attrs
import narwhals
import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    FLIP_NEGATE,
    ExtraColumnRule,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_HAZARD_RATIO_95L_COL,
    GWASLAB_HAZARD_RATIO_95U_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_ODDS_RATIO_95L_COL,
    GWASLAB_ODDS_RATIO_95U_COL,
    GWASLAB_ODDS_RATIO_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.regenie_constants import REGENIE_A1FREQ_CASES_COL
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    CONSISTENT_INDEL,
    CONSISTENT_SNV,
    INCONSISTENT_SNV,
    TEST_OPTIONS,
    Variant,
    positions,
    row_at,
    run_harmonization,
    sumstats_frame,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
EAF = GWASLAB_EFFECT_ALLELE_FREQ_COL
BETA = GWASLAB_BETA_COL
OR_95L = GWASLAB_ODDS_RATIO_95L_COL
OR_95U = GWASLAB_ODDS_RATIO_95U_COL
HR_95L = GWASLAB_HAZARD_RATIO_95L_COL
HR_95U = GWASLAB_HAZARD_RATIO_95U_COL
_UNREGISTERED_COLUMN = "MYSTERY_STATISTIC"
_EXTRA_COLUMN = "BETA_ALTERNATIVE_MODEL"


def test_snvs_are_oriented_against_the_reference(tmp_path: Path) -> None:
    variants = [
        CONSISTENT_SNV,  # 1: NEA is the reference base -> kept
        INCONSISTENT_SNV,  # 2: EA is the reference base -> swapped
        Variant(pos=3, ea="A", nea="C", beta=0.3),  # complement of NEA is ref -> complemented
        Variant(pos=11, ea="G", nea="A", beta=0.4),  # complement of EA is ref -> complemented and swapped
        Variant(pos=12, ea="C", nea="G"),  # nothing matches -> dropped
    ]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert result.select(GWASLAB_POS_COL, EA, NEA, BETA, EAF).rows() == [
        (1, "G", "A", pytest.approx(0.1), pytest.approx(0.3)),
        (2, "T", "C", pytest.approx(-0.2), pytest.approx(0.7)),
        (3, "T", "G", pytest.approx(0.3), pytest.approx(0.3)),
        (11, "T", "C", pytest.approx(-0.4), pytest.approx(0.7)),
    ]


def test_indels_with_one_matching_allele_are_oriented_by_it(tmp_path: Path) -> None:
    variants = [
        CONSISTENT_SNV,
        CONSISTENT_INDEL,  # 21: only NEA G matches -> kept
        Variant(pos=22, ea="A", nea="AG", beta=0.5),  # only EA matches -> swapped
        Variant(pos=23, ea="C", nea="CA"),  # neither matches -> dropped
    ]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert positions(result) == [1, 21, 22]
    swapped = row_at(result, 22)
    assert (swapped[EA], swapped[NEA]) == ("AG", "A")
    assert swapped[BETA] == pytest.approx(-0.5)


def test_trusted_table_keeps_ambiguous_indel_in_source_orientation(tmp_path: Path) -> None:
    ambiguous = Variant(pos=5, ea="T", nea="TG")  # both alleles match the T-G-run
    result = run_harmonization(
        tmp_path / "run", sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, ambiguous])
    )
    kept = row_at(result, 5)
    assert (kept[EA], kept[NEA]) == ("T", "TG")


@pytest.mark.parametrize("trusted", [True, False])
def test_palindromic_mnp_is_kept_only_when_trusted(tmp_path: Path, trusted: bool) -> None:
    palindromic_mnp = Variant(pos=41, ea="GT", nea="AC")
    variants = [CONSISTENT_SNV, CONSISTENT_INDEL, palindromic_mnp]
    if not trusted:
        variants.append(INCONSISTENT_SNV)
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert (41 in positions(result)) == trusted


def test_invalid_alleles_are_dropped_and_lowercase_is_accepted(tmp_path: Path) -> None:
    variants = [
        CONSISTENT_SNV,
        Variant(pos=3, ea="N", nea="C"),
        Variant(pos=4, ea="T", nea="T"),
        Variant(pos=51, ea="a", nea="g"),
    ]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants))
    assert positions(result) == [1, 51]
    assert (row_at(result, 51)[EA], row_at(result, 51)[NEA]) == ("A", "G")


def test_flip_covers_allele_frequency_columns_and_confidence_bounds(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV]).with_columns(
        pl.Series(REGENIE_A1FREQ_CASES_COL, [0.3, 0.25]),
        pl.Series(GWASLAB_ODDS_RATIO_COL, [2.0, 4.0]),
        pl.Series(OR_95L, [1.5, 2.0]),
        pl.Series(OR_95U, [2.5, 8.0]),
        pl.Series(HR_95L, [1.2, 1.25]),
        pl.Series(HR_95U, [1.8, 5.0]),
    )
    result = run_harmonization(tmp_path / "run", frame)
    unchanged = row_at(result, 1)
    assert (unchanged[OR_95L], unchanged[OR_95U]) == (pytest.approx(1.5), pytest.approx(2.5))
    swapped = row_at(result, 2)
    assert swapped[REGENIE_A1FREQ_CASES_COL] == pytest.approx(0.75)
    assert swapped[GWASLAB_ODDS_RATIO_COL] == pytest.approx(0.25)
    assert (swapped[OR_95L], swapped[OR_95U]) == (pytest.approx(0.125), pytest.approx(0.5))
    assert (swapped[HR_95L], swapped[HR_95U]) == (pytest.approx(0.2), pytest.approx(0.8))


def test_unregistered_column_fails(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV]).with_columns(pl.lit(1.0).alias(_UNREGISTERED_COLUMN))
    with pytest.raises(AssertionError):
        run_harmonization(tmp_path / "run", frame)


def test_extra_column_rule_is_applied(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV]).with_columns(
        pl.Series(_EXTRA_COLUMN, [0.5, 0.5])
    )
    options = attrs.evolve(
        TEST_OPTIONS,
        extra_column_rules=(ExtraColumnRule(column=_EXTRA_COLUMN, rule=FLIP_NEGATE),),
    )
    result = run_harmonization(tmp_path / "run", frame, options=options)
    assert row_at(result, 1)[_EXTRA_COLUMN] == pytest.approx(0.5)
    assert row_at(result, 2)[_EXTRA_COLUMN] == pytest.approx(-0.5)


def test_output_does_not_depend_on_input_row_order(tmp_path: Path) -> None:
    chr1_keep, chr1_swap = CONSISTENT_SNV, INCONSISTENT_SNV
    chr2_keep = Variant(chrom=2, pos=5, ea="G", nea="A")
    chr2_swap = Variant(chrom=2, pos=9, ea="C", nea="T")
    ordered = run_harmonization(
        tmp_path / "ordered", sumstats_frame([chr1_keep, chr1_swap, chr2_keep, chr2_swap])
    )
    interleaved = run_harmonization(
        tmp_path / "interleaved",
        sumstats_frame([chr2_swap, chr1_keep, chr2_keep, chr1_swap]),
    )
    assert interleaved.equals(ordered)
    assert ordered[GWASLAB_CHROM_COL].to_list() == [1, 1, 2, 2]


def test_excluded_chromosome_rows_are_dropped(tmp_path: Path) -> None:
    mitochondrial = Variant(chrom=25, pos=1, ea="G", nea="A")
    result = run_harmonization(tmp_path / "run", sumstats_frame([CONSISTENT_SNV, mitochondrial]))
    assert result[GWASLAB_CHROM_COL].to_list() == [1]


def test_chromosome_missing_from_fasta_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        run_harmonization(
            tmp_path / "run",
            sumstats_frame([CONSISTENT_SNV, Variant(chrom=3, pos=1, ea="G", nea="A")]),
        )


def test_null_position_fails(tmp_path: Path) -> None:
    frame = sumstats_frame([CONSISTENT_SNV]).with_columns(
        pl.lit(None, dtype=pl.Int64).alias(GWASLAB_POS_COL)
    )
    with pytest.raises(AssertionError):
        run_harmonization(tmp_path / "run", frame)


def test_duplicate_variant_after_orientation_fails(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, Variant(pos=2, ea="T", nea="C")]
    with pytest.raises(AssertionError):
        run_harmonization(tmp_path / "run", sumstats_frame(variants))


def test_long_alleles_are_classified_with_a_tiny_gather_budget(tmp_path: Path) -> None:
    long_mnp = Variant(pos=31, ea="T" * 10, nea="T" * 9 + "A", beta=0.2)  # EA matches -> swap
    options = attrs.evolve(TEST_OPTIONS, max_gather_bytes=1)
    result = run_harmonization(
        tmp_path / "run", sumstats_frame([CONSISTENT_SNV, long_mnp]), options=options
    )
    swapped = row_at(result, 31)
    assert (swapped[EA], swapped[NEA]) == ("T" * 9 + "A", "T" * 10)
    assert swapped[BETA] == pytest.approx(-0.2)


@frozen
class _SwitchToDuckDbPipe(DataProcessingPipe):
    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.collect().lazy(backend="duckdb")


def test_pipe_that_changes_backend_fails(tmp_path: Path) -> None:
    with pytest.raises(AssertionError):
        run_harmonization(
            tmp_path / "run", sumstats_frame([CONSISTENT_SNV]), pipe=_SwitchToDuckDbPipe()
        )
```

- [ ] **Step 3: Run the tests to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py -q`
Expected: collection error, ModuleNotFoundError for flip or genome_reference_harmonization_task.

- [ ] **Step 4: Add column-name constants**

Append to mecfs_bio/constants/gwaslab_constants.py, below GWASLAB_SNPID_COL:

```python
# Further gwaslab standard column names, from the "gwaslab" entry of
# https://github.com/Cloufield/formatbook/blob/main/formatbook.json
GWASLAB_Z_COL = "Z"
GWASLAB_T_STATISTIC_COL = "T"
GWASLAB_F_STATISTIC_COL = "F"
GWASLAB_HAZARD_RATIO_COL = "HR"
GWASLAB_ODDS_RATIO_95L_COL = "OR_95L"
GWASLAB_ODDS_RATIO_95U_COL = "OR_95U"
GWASLAB_HAZARD_RATIO_95L_COL = "HR_95L"
GWASLAB_HAZARD_RATIO_95U_COL = "HR_95U"
GWASLAB_P_HET_COL = "P_HET"
GWASLAB_I2_COL = "I2"
GWASLAB_SNPR2_COL = "SNPR2"
GWASLAB_DOF_COL = "DOF"
GWASLAB_MAF_COL = "MAF"
```

Append to mecfs_bio/constants/regenie_constants.py:

```python
# Binary-trait output columns. gwaslab passes them through without renaming, so they
# survive into gwaslab-format tables (DecodeME carries all of them).
REGENIE_A1FREQ_CASES_COL = "A1FREQ_CASES"  # frequency of ALLELE1 in cases
REGENIE_A1FREQ_CONTROLS_COL = "A1FREQ_CONTROLS"  # frequency of ALLELE1 in controls
REGENIE_N_CASES_COL = "N_CASES"
REGENIE_N_CONTROLS_COL = "N_CONTROLS"
REGENIE_TEST_COL = "TEST"
REGENIE_EXTRA_COL = "EXTRA"
```

- [ ] **Step 5: Implement flip.py**

Registry scope:
- **Included:** gwaslab's standard column names (the formatbook "gwaslab" entry) plus the columns this repo's datasets actually carry: the regenie binary-trait columns, and N_EFF from GWASLAB_EFFECTIVE_SAMPLE_SIZE.
- **Not registered:** DIRECTION (a meta-analysis direction string). Only formats this repo never uses map to it (metal, mrmega, the auto formats), so it fails loudly if it ever appears.
- **BETA_95L/BETA_95U** are not gwaslab standard columns and are not used here.
- **REF and ALT** are standard, but not registered (decided in plan review, 2026-09-15). After harmonization NEA and EA carry the reference orientation, and a source REF/ALT pair could contradict it, so an input carrying them fails loudly.

```python
"""
How each summary-statistic column changes when a variant's effect allele is swapped.

Columns are expected to use gwaslab's standard names (the "gwaslab" entry of formatbook.json), plus the regenie binary-trait columns and effective sample size used by datasets in this repo. Every non-allele column must be registered here, or declared through ExtraColumnRule, so that a new allele-dependent statistic can never be silently left unflipped.

All flipped values are computed from the original frame in a single with_columns call, so a confidence bound is never read after being overwritten.

FlipRule is what the type checker sees and the FLIP_* constants are what code uses. An import-time assertion keeps the two in step. A StrEnum cannot replace them, because polars converts Python Enum members into its own Enum dtype.
"""

from collections.abc import Mapping, Sequence
from typing import Final, Literal, get_args

import polars as pl
from attrs import frozen

from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHISQ_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_DOF_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_EFFECTIVE_SAMPLE_SIZE,
    GWASLAB_F_STATISTIC_COL,
    GWASLAB_HAZARD_RATIO_95L_COL,
    GWASLAB_HAZARD_RATIO_95U_COL,
    GWASLAB_HAZARD_RATIO_COL,
    GWASLAB_I2_COL,
    GWASLAB_INFO_SCORE_COL,
    GWASLAB_MAF_COL,
    GWASLAB_MLOG10P_COL,
    GWASLAB_N_CASE_COL,
    GWASLAB_N_CONTROL_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_ODDS_RATIO_95L_COL,
    GWASLAB_ODDS_RATIO_95U_COL,
    GWASLAB_ODDS_RATIO_COL,
    GWASLAB_P_COL,
    GWASLAB_P_HET_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
    GWASLAB_SNPID_COL,
    GWASLAB_SNPR2_COL,
    GWASLAB_STATUS_COL,
    GWASLAB_T_STATISTIC_COL,
    GWASLAB_Z_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_CASES_COL,
    REGENIE_A1FREQ_CONTROLS_COL,
    REGENIE_EXTRA_COL,
    REGENIE_N_CASES_COL,
    REGENIE_N_CONTROLS_COL,
    REGENIE_TEST_COL,
)

FlipRule = Literal["negate", "complement", "invert", "invariant"]
FLIP_NEGATE: Final[FlipRule] = "negate"
FLIP_COMPLEMENT: Final[FlipRule] = "complement"
FLIP_INVERT: Final[FlipRule] = "invert"
FLIP_INVARIANT: Final[FlipRule] = "invariant"
assert set(get_args(FlipRule)) == {FLIP_NEGATE, FLIP_COMPLEMENT, FLIP_INVERT, FLIP_INVARIANT}


@frozen
class InvertedBoundPair:
    """Ratio confidence bounds: on a flip, lower becomes 1/upper and upper becomes 1/lower."""

    lower: str
    upper: str


@frozen
class ExtraColumnRule:
    column: str
    rule: FlipRule


COLUMN_FLIP_RULES: Mapping[str, FlipRule] = {
    GWASLAB_BETA_COL: FLIP_NEGATE,
    GWASLAB_Z_COL: FLIP_NEGATE,
    GWASLAB_T_STATISTIC_COL: FLIP_NEGATE,
    GWASLAB_EFFECT_ALLELE_FREQ_COL: FLIP_COMPLEMENT,
    REGENIE_A1FREQ_CASES_COL: FLIP_COMPLEMENT,
    REGENIE_A1FREQ_CONTROLS_COL: FLIP_COMPLEMENT,
    GWASLAB_ODDS_RATIO_COL: FLIP_INVERT,
    GWASLAB_HAZARD_RATIO_COL: FLIP_INVERT,
    GWASLAB_SNPID_COL: FLIP_INVARIANT,
    GWASLAB_RSID_COL: FLIP_INVARIANT,
    GWASLAB_CHROM_COL: FLIP_INVARIANT,
    GWASLAB_POS_COL: FLIP_INVARIANT,
    GWASLAB_SE_COL: FLIP_INVARIANT,
    GWASLAB_P_COL: FLIP_INVARIANT,
    GWASLAB_MLOG10P_COL: FLIP_INVARIANT,
    GWASLAB_CHISQ_COL: FLIP_INVARIANT,
    GWASLAB_F_STATISTIC_COL: FLIP_INVARIANT,
    GWASLAB_P_HET_COL: FLIP_INVARIANT,
    GWASLAB_I2_COL: FLIP_INVARIANT,
    GWASLAB_SNPR2_COL: FLIP_INVARIANT,
    GWASLAB_DOF_COL: FLIP_INVARIANT,
    GWASLAB_SAMPLE_SIZE_COLUMN: FLIP_INVARIANT,
    GWASLAB_N_CASE_COL: FLIP_INVARIANT,
    GWASLAB_N_CONTROL_COL: FLIP_INVARIANT,
    REGENIE_N_CASES_COL: FLIP_INVARIANT,
    REGENIE_N_CONTROLS_COL: FLIP_INVARIANT,
    GWASLAB_EFFECTIVE_SAMPLE_SIZE: FLIP_INVARIANT,
    GWASLAB_INFO_SCORE_COL: FLIP_INVARIANT,
    GWASLAB_MAF_COL: FLIP_INVARIANT,
    REGENIE_TEST_COL: FLIP_INVARIANT,
    REGENIE_EXTRA_COL: FLIP_INVARIANT,
}

BOUND_PAIRS: tuple[InvertedBoundPair, ...] = (
    InvertedBoundPair(lower=GWASLAB_ODDS_RATIO_95L_COL, upper=GWASLAB_ODDS_RATIO_95U_COL),
    InvertedBoundPair(lower=GWASLAB_HAZARD_RATIO_95L_COL, upper=GWASLAB_HAZARD_RATIO_95U_COL),
)

# Dropped from the output: STATUS describes gwaslab's processing, not this Task's.
DROPPED_COLUMNS: frozenset[str] = frozenset({GWASLAB_STATUS_COL})
_ALLELE_COLUMNS = frozenset({GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL})


def resolve_column_rules(
    columns: Sequence[str], extra: Sequence[ExtraColumnRule]
) -> Mapping[str, FlipRule]:
    """Flip rule for every present non-allele column; asserts that none is unregistered."""
    extra_rules = {rule.column: rule.rule for rule in extra}
    assert len(extra_rules) == len(extra), "extra_column_rules names a column twice"
    rules = {**COLUMN_FLIP_RULES, **extra_rules}
    bound_columns = {name for pair in BOUND_PAIRS for name in (pair.lower, pair.upper)}
    unregistered = [
        column
        for column in columns
        if column not in rules
        and column not in bound_columns
        and column not in _ALLELE_COLUMNS
        and column not in DROPPED_COLUMNS
    ]
    assert not unregistered, (
        f"columns without a flip rule: {unregistered}; register them in flip.py or pass "
        "ExtraColumnRule entries"
    )
    for pair in BOUND_PAIRS:
        assert (pair.lower in columns) == (pair.upper in columns), (
            f"confidence bounds {pair.lower} and {pair.upper} must be present together"
        )
    return {column: rules[column] for column in columns if column in rules}


def _flipped(column: str, rule: FlipRule) -> pl.Expr:
    value = pl.col(column)
    if rule == FLIP_NEGATE:
        return -value
    if rule == FLIP_COMPLEMENT:
        return 1 - value
    assert rule == FLIP_INVERT, f"rule {rule} has no flipped expression"
    return 1 / value


def flip_statistics(
    frame: pl.DataFrame, mask: pl.Expr, rules: Mapping[str, FlipRule]
) -> pl.DataFrame:
    """Apply every flip rule to the rows selected by mask; alleles are not touched."""
    updates = [
        pl.when(mask).then(_flipped(column, rule)).otherwise(pl.col(column)).alias(column)
        for column, rule in rules.items()
        if rule != FLIP_INVARIANT and column in frame.columns
    ]
    for pair in BOUND_PAIRS:
        if pair.lower in frame.columns:
            updates.append(
                pl.when(mask).then(1 / pl.col(pair.upper)).otherwise(pl.col(pair.lower)).alias(pair.lower)
            )
            updates.append(
                pl.when(mask).then(1 / pl.col(pair.lower)).otherwise(pl.col(pair.upper)).alias(pair.upper)
            )
    return frame.with_columns(updates) if updates else frame
```

- [ ] **Step 6: Implement options.py and outcomes.py**

options.py:

```python
"""Options for genome-reference harmonization, validated at construction."""

from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    DEFAULT_MAX_GATHER_BYTES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    ExtraColumnRule,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_CODE_FOR_NAME


@frozen
class GenomeReferenceHarmonizationOptions:
    """
    min_checkable_snvs, min_checkable_indels: a table is trusted only if it has at least
        this many checkable SNVs (non-palindromic) and indels, all consistent with NEA as the
        reference allele.
    palindrome_maf_threshold, panel_maf_threshold: an untrusted palindromic SNV is resolved
        by frequency only when both its summary-statistic MAF and its panel MAF are at most
        these (gwaslab's defaults).
    indel_max_af_distance, indel_min_af_margin: stringent rules for untrusted ambiguous indels.
        A reading is chosen only if its predicted EAF is within the distance and beats the
        other reading by at least the margin.
    keep_unresolved_palindromes: keep untrusted palindromic SNVs whose strand cannot be
        resolved, instead of dropping them. Never applies to indels.
    excluded_chromosomes: gwaslab chromosome codes whose rows are dropped (MT by default,
        since UCSC hg19 chrM is not rCRS).
    extra_column_rules: flip rules for dataset-specific columns.
    max_gather_bytes: memory budget for one reference-gather slice.
    """

    min_checkable_snvs: int = 10_000
    min_checkable_indels: int = 1_000
    palindrome_maf_threshold: float = 0.4
    panel_maf_threshold: float = 0.4
    indel_max_af_distance: float = 0.1
    indel_min_af_margin: float = 0.2
    keep_unresolved_palindromes: bool = False
    excluded_chromosomes: tuple[int, ...] = (GWASLAB_CHROM_CODE_FOR_NAME["MT"],)
    extra_column_rules: tuple[ExtraColumnRule, ...] = ()
    max_gather_bytes: int = DEFAULT_MAX_GATHER_BYTES

    def __attrs_post_init__(self) -> None:
        assert self.min_checkable_snvs >= 1 and self.min_checkable_indels >= 1
        assert 0 < self.palindrome_maf_threshold < 0.5
        assert 0 < self.panel_maf_threshold < 0.5
        assert 0 < self.indel_max_af_distance < 1
        assert 0 < self.indel_min_af_margin < 1, (
            "a strictly positive margin keeps the keep and flip readings from both being chosen"
        )
        assert self.max_gather_bytes > 0
```

outcomes.py:

```python
"""
Named outcomes of genome-reference harmonization: allele actions, palindrome decisions
and drop reasons.

Each Literal alias is what the type checker sees; the constants are what code uses, so
no outcome string is typed more than twice (once in the alias, once in its constant).
An import-time assertion keeps each alias and its constants in step. A StrEnum cannot
replace them, because polars converts Python Enum members into its own Enum dtype
inside expressions.
"""

from typing import Final, Literal, get_args

AlleleAction = Literal["keep", "swap", "complement", "complement_swap"]
ACTION_KEEP: Final[AlleleAction] = "keep"
ACTION_SWAP: Final[AlleleAction] = "swap"
ACTION_COMPLEMENT: Final[AlleleAction] = "complement"
ACTION_COMPLEMENT_SWAP: Final[AlleleAction] = "complement_swap"
assert set(get_args(AlleleAction)) == {
    ACTION_KEEP,
    ACTION_SWAP,
    ACTION_COMPLEMENT,
    ACTION_COMPLEMENT_SWAP,
}

PalindromeDecision = Literal["keep", "strand_flip", "unresolved"]
PALINDROME_KEEP: Final[PalindromeDecision] = "keep"
PALINDROME_STRAND_FLIP: Final[PalindromeDecision] = "strand_flip"
PALINDROME_UNRESOLVED: Final[PalindromeDecision] = "unresolved"
assert set(get_args(PalindromeDecision)) == {
    PALINDROME_KEEP,
    PALINDROME_STRAND_FLIP,
    PALINDROME_UNRESOLVED,
}

DropReason = Literal[
    "invalid_allele",
    "not_on_reference",
    "indel_not_on_reference",
    "palindromic_mnp_untrusted",
    "palindrome_unresolved",
    "ambiguous_indel_no_eaf",
    "ambiguous_indel_not_in_panel",
    "ambiguous_indel_af_mismatch",
    "ambiguous_indel_af_indecisive",
]
DROP_INVALID_ALLELE: Final[DropReason] = "invalid_allele"
DROP_NOT_ON_REFERENCE: Final[DropReason] = "not_on_reference"
DROP_INDEL_NOT_ON_REFERENCE: Final[DropReason] = "indel_not_on_reference"
DROP_PALINDROMIC_MNP_UNTRUSTED: Final[DropReason] = "palindromic_mnp_untrusted"
DROP_PALINDROME_UNRESOLVED: Final[DropReason] = "palindrome_unresolved"
DROP_AMBIGUOUS_INDEL_NO_EAF: Final[DropReason] = "ambiguous_indel_no_eaf"
DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL: Final[DropReason] = "ambiguous_indel_not_in_panel"
DROP_AMBIGUOUS_INDEL_AF_MISMATCH: Final[DropReason] = "ambiguous_indel_af_mismatch"
DROP_AMBIGUOUS_INDEL_AF_INDECISIVE: Final[DropReason] = "ambiguous_indel_af_indecisive"
assert set(get_args(DropReason)) == {
    DROP_INVALID_ALLELE,
    DROP_NOT_ON_REFERENCE,
    DROP_INDEL_NOT_ON_REFERENCE,
    DROP_PALINDROMIC_MNP_UNTRUSTED,
    DROP_PALINDROME_UNRESOLVED,
    DROP_AMBIGUOUS_INDEL_NO_EAF,
    DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL,
    DROP_AMBIGUOUS_INDEL_AF_MISMATCH,
    DROP_AMBIGUOUS_INDEL_AF_INDECISIVE,
}
```

- [ ] **Step 7: Implement allele_classes.py**

```python
"""
Per-row allele classes relative to the genome reference.

Equal-length variants (SNVs, MNPs) are compared on the plus strand. The reverse
complement is consulted only when neither plus-strand allele matches, as gwaslab does.
Different-length variants (indels) are compared on the plus strand only: a
reverse-complemented left-anchored indel loses its anchor base, so a minus-strand match
would be spurious.

AlleleClass is what the type checker sees; the CLASS_* constants are what code uses. An
import-time assertion keeps the two in step (see outcomes.py for why this is not a
StrEnum).
"""

from typing import Final, Literal, get_args

import numpy as np
import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

AlleleClass = Literal[
    "nea_ref",
    "ea_ref",
    "nea_ref_rc",
    "ea_ref_rc",
    "not_on_reference",
    "indel_both",
    "indel_nea_only",
    "indel_ea_only",
    "indel_not_on_reference",
]
CLASS_NEA_REF: Final[AlleleClass] = "nea_ref"
CLASS_EA_REF: Final[AlleleClass] = "ea_ref"
CLASS_NEA_REF_RC: Final[AlleleClass] = "nea_ref_rc"
CLASS_EA_REF_RC: Final[AlleleClass] = "ea_ref_rc"
CLASS_NOT_ON_REFERENCE: Final[AlleleClass] = "not_on_reference"
CLASS_INDEL_BOTH: Final[AlleleClass] = "indel_both"
CLASS_INDEL_NEA_ONLY: Final[AlleleClass] = "indel_nea_only"
CLASS_INDEL_EA_ONLY: Final[AlleleClass] = "indel_ea_only"
CLASS_INDEL_NOT_ON_REFERENCE: Final[AlleleClass] = "indel_not_on_reference"
assert set(get_args(AlleleClass)) == {
    CLASS_NEA_REF,
    CLASS_EA_REF,
    CLASS_NEA_REF_RC,
    CLASS_EA_REF_RC,
    CLASS_NOT_ON_REFERENCE,
    CLASS_INDEL_BOTH,
    CLASS_INDEL_NEA_ONLY,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
}

ALLELE_CLASS_COL = "_allele_class"
IS_PALINDROMIC_SNV_COL = "_is_palindromic_snv"
IS_PALINDROMIC_MNP_COL = "_is_palindromic_mnp"
_VALID_ALLELE_PATTERN = "^[ACGT]+$"
_NEA_RC_COL = "_nea_rc"
_EA_RC_COL = "_ea_rc"
_BASES = ["A", "C", "G", "T"]
_COMPLEMENTS = ["T", "G", "C", "A"]


def reverse_complement_expr(column: str) -> pl.Expr:
    return pl.col(column).str.reverse().str.replace_many(_BASES, _COMPLEMENTS)


def prepare_alleles(frame: pl.DataFrame) -> pl.DataFrame:
    """Alleles as uppercase strings and positions as Int64."""
    return frame.with_columns(
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String).str.to_uppercase(),
        pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).cast(pl.String).str.to_uppercase(),
        pl.col(GWASLAB_POS_COL).cast(pl.Int64),
    )


def valid_alleles_expr() -> pl.Expr:
    ea = pl.col(GWASLAB_EFFECT_ALLELE_COL)
    nea = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL)
    return (
        ea.str.contains(_VALID_ALLELE_PATTERN)
        & nea.str.contains(_VALID_ALLELE_PATTERN)
        & (ea != nea)
    )


def classify_alleles(
    frame: pl.DataFrame, fasta: IndexedFasta, chrom: int, max_gather_bytes: int
) -> pl.DataFrame:
    """Add ALLELE_CLASS_COL and the palindrome flags to rows with valid, uppercase alleles."""
    pos = frame[GWASLAB_POS_COL].to_numpy()
    nea_match = reference_matches(
        fasta,
        chrom=chrom,
        positions=pos,
        alleles=frame[GWASLAB_NON_EFFECT_ALLELE_COL],
        max_gather_bytes=max_gather_bytes,
    )
    ea_match = reference_matches(
        fasta,
        chrom=chrom,
        positions=pos,
        alleles=frame[GWASLAB_EFFECT_ALLELE_COL],
        max_gather_bytes=max_gather_bytes,
    )
    equal_length = (
        frame[GWASLAB_EFFECT_ALLELE_COL].str.len_bytes()
        == frame[GWASLAB_NON_EFFECT_ALLELE_COL].str.len_bytes()
    ).to_numpy()
    needs_rc = equal_length & ~nea_match & ~ea_match
    rc_nea_match = np.zeros(frame.height, dtype=bool)
    rc_ea_match = np.zeros(frame.height, dtype=bool)
    rc_rows = np.flatnonzero(needs_rc)
    if len(rc_rows):
        rc = frame.select(
            reverse_complement_expr(GWASLAB_NON_EFFECT_ALLELE_COL).alias(_NEA_RC_COL),
            reverse_complement_expr(GWASLAB_EFFECT_ALLELE_COL).alias(_EA_RC_COL),
        )
        rc_nea_match[rc_rows] = reference_matches(
            fasta,
            chrom=chrom,
            positions=pos[rc_rows],
            alleles=rc[_NEA_RC_COL].gather(rc_rows),
            max_gather_bytes=max_gather_bytes,
        )
        rc_ea_match[rc_rows] = reference_matches(
            fasta,
            chrom=chrom,
            positions=pos[rc_rows],
            alleles=rc[_EA_RC_COL].gather(rc_rows),
            max_gather_bytes=max_gather_bytes,
        )
    classes = np.select(
        [
            equal_length & nea_match,
            equal_length & ea_match,
            needs_rc & rc_nea_match,
            needs_rc & rc_ea_match,
            equal_length,
            nea_match & ea_match,
            nea_match,
            ea_match,
        ],
        [
            CLASS_NEA_REF,
            CLASS_EA_REF,
            CLASS_NEA_REF_RC,
            CLASS_EA_REF_RC,
            CLASS_NOT_ON_REFERENCE,
            CLASS_INDEL_BOTH,
            CLASS_INDEL_NEA_ONLY,
            CLASS_INDEL_EA_ONLY,
        ],
        default=CLASS_INDEL_NOT_ON_REFERENCE,
    )
    ea_length = pl.col(GWASLAB_EFFECT_ALLELE_COL).str.len_bytes()
    same_length = ea_length == pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).str.len_bytes()
    palindromic = same_length & (
        reverse_complement_expr(GWASLAB_NON_EFFECT_ALLELE_COL)
        == pl.col(GWASLAB_EFFECT_ALLELE_COL)
    )
    return frame.with_columns(
        pl.Series(ALLELE_CLASS_COL, classes, dtype=pl.String),
        (palindromic & (ea_length == 1)).alias(IS_PALINDROMIC_SNV_COL),
        (palindromic & (ea_length > 1)).alias(IS_PALINDROMIC_MNP_COL),
    )
```

- [ ] **Step 8: Implement trust.py**

```python
"""
Whether a table's NEA column is trustworthy as the reference allele.

Trust requires exactly 100% consistency. That is, every checkable SNV (non-palindromic,
single base) has NEA as the reference base, and every checkable indel (different
lengths) has only NEA matching the reference. Each set must also meet a minimum count.
Ambiguous indels, whose alleles both match, carry no evidence and are not counted.
"""

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NEA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
    CLASS_NEA_REF,
    IS_PALINDROMIC_SNV_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
)

_CONSISTENT_SNVS = "consistent_snvs"
_INCONSISTENT_SNVS = "inconsistent_snvs"
_CONSISTENT_INDELS = "consistent_indels"
_INCONSISTENT_INDELS = "inconsistent_indels"


@frozen
class TrustCounts:
    consistent_snvs: int
    inconsistent_snvs: int
    consistent_indels: int
    inconsistent_indels: int

    def __add__(self, other: "TrustCounts") -> "TrustCounts":
        return TrustCounts(
            consistent_snvs=self.consistent_snvs + other.consistent_snvs,
            inconsistent_snvs=self.inconsistent_snvs + other.inconsistent_snvs,
            consistent_indels=self.consistent_indels + other.consistent_indels,
            inconsistent_indels=self.inconsistent_indels + other.inconsistent_indels,
        )

    @classmethod
    def zero(cls) -> "TrustCounts":
        return cls(
            consistent_snvs=0, inconsistent_snvs=0, consistent_indels=0, inconsistent_indels=0
        )


def count_trust_evidence(classified: pl.DataFrame) -> TrustCounts:
    """Count checkable SNVs and indels in a frame produced by classify_alleles."""
    allele_class = pl.col(ALLELE_CLASS_COL)
    ea_length = pl.col(GWASLAB_EFFECT_ALLELE_COL).str.len_bytes()
    nea_length = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL).str.len_bytes()
    checkable_snv = (ea_length == 1) & (nea_length == 1) & ~pl.col(IS_PALINDROMIC_SNV_COL)
    indel = ea_length != nea_length
    counts = classified.select(
        (checkable_snv & (allele_class == CLASS_NEA_REF)).sum().alias(_CONSISTENT_SNVS),
        (checkable_snv & (allele_class != CLASS_NEA_REF)).sum().alias(_INCONSISTENT_SNVS),
        (indel & (allele_class == CLASS_INDEL_NEA_ONLY)).sum().alias(_CONSISTENT_INDELS),
        (indel & allele_class.is_in([CLASS_INDEL_EA_ONLY, CLASS_INDEL_NOT_ON_REFERENCE]))
        .sum()
        .alias(_INCONSISTENT_INDELS),
    ).row(0, named=True)
    return TrustCounts(
        consistent_snvs=int(counts[_CONSISTENT_SNVS]),
        inconsistent_snvs=int(counts[_INCONSISTENT_SNVS]),
        consistent_indels=int(counts[_CONSISTENT_INDELS]),
        inconsistent_indels=int(counts[_INCONSISTENT_INDELS]),
    )


def decide_trust(counts: TrustCounts, options: GenomeReferenceHarmonizationOptions) -> bool:
    return (
        counts.inconsistent_snvs == 0
        and counts.inconsistent_indels == 0
        and counts.consistent_snvs >= options.min_checkable_snvs
        and counts.consistent_indels >= options.min_checkable_indels
    )
```

- [ ] **Step 9: Implement resolve_chromosome.py**

```python
"""
Resolve one chromosome of summary statistics against the genome reference.

Every input row is returned. Alleles are oriented so that NEA is the plus-strand
reference allele, allele-dependent statistics are flipped to match, and
DROP_REASON_COL is null for rows to keep. The function holds one chromosome in memory.
"""

from collections.abc import Mapping
from typing import Protocol

import polars as pl
from attrs import frozen

from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_EA_REF,
    CLASS_EA_REF_RC,
    CLASS_INDEL_EA_ONLY,
    CLASS_INDEL_NOT_ON_REFERENCE,
    CLASS_NEA_REF_RC,
    CLASS_NOT_ON_REFERENCE,
    IS_PALINDROMIC_MNP_COL,
    classify_alleles,
    prepare_alleles,
    reverse_complement_expr,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    FlipRule,
    flip_statistics,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_COMPLEMENT,
    ACTION_COMPLEMENT_SWAP,
    ACTION_KEEP,
    ACTION_SWAP,
    DROP_INDEL_NOT_ON_REFERENCE,
    DROP_INVALID_ALLELE,
    DROP_NOT_ON_REFERENCE,
    DROP_PALINDROMIC_MNP_UNTRUSTED,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
)

DROP_REASON_COL = "_drop_reason"
ROW_INDEX_COL = "_row_index"
ACTION_COL = "_action"


class PanelLoader(Protocol):
    def __call__(self, positions: pl.Series) -> pl.DataFrame:
        """Panel rows (POS Int64, REF, ALT, AF) of the current chromosome at these positions."""
        ...


@frozen
class ChromosomeContext:
    chrom: int
    fasta: IndexedFasta
    trusted: bool
    rules: Mapping[str, FlipRule]
    options: GenomeReferenceHarmonizationOptions


def resolve_chromosome(
    rows: pl.DataFrame, context: ChromosomeContext, load_panel: PanelLoader
) -> pl.DataFrame:
    output_columns = [*rows.columns, DROP_REASON_COL]
    frame = prepare_alleles(rows).with_row_index(ROW_INDEX_COL)
    invalid = frame.filter(~valid_alleles_expr()).with_columns(
        pl.lit(DROP_INVALID_ALLELE, dtype=pl.String).alias(DROP_REASON_COL)
    )
    valid = classify_alleles(
        frame.filter(valid_alleles_expr()),
        fasta=context.fasta,
        chrom=context.chrom,
        max_gather_bytes=context.options.max_gather_bytes,
    ).with_columns(_base_action_expr(), _base_drop_reason_expr(trusted=context.trusted))
    valid = _apply_allele_actions(valid, context.rules)
    return (
        pl.concat(
            [
                valid.select(ROW_INDEX_COL, *output_columns),
                invalid.select(ROW_INDEX_COL, *output_columns),
            ]
        )
        .sort(ROW_INDEX_COL)
        .drop(ROW_INDEX_COL)
    )


def _base_action_expr() -> pl.Expr:
    allele_class = pl.col(ALLELE_CLASS_COL)
    return (
        pl.when(allele_class.is_in([CLASS_EA_REF, CLASS_INDEL_EA_ONLY]))
        .then(pl.lit(ACTION_SWAP))
        .when(allele_class == CLASS_NEA_REF_RC)
        .then(pl.lit(ACTION_COMPLEMENT))
        .when(allele_class == CLASS_EA_REF_RC)
        .then(pl.lit(ACTION_COMPLEMENT_SWAP))
        .otherwise(pl.lit(ACTION_KEEP))
        .alias(ACTION_COL)
    )


def _base_drop_reason_expr(trusted: bool) -> pl.Expr:
    allele_class = pl.col(ALLELE_CLASS_COL)
    untrusted_palindromic_mnp = pl.col(IS_PALINDROMIC_MNP_COL) & pl.lit(not trusted)
    return (
        pl.when(allele_class == CLASS_NOT_ON_REFERENCE)
        .then(pl.lit(DROP_NOT_ON_REFERENCE))
        .when(allele_class == CLASS_INDEL_NOT_ON_REFERENCE)
        .then(pl.lit(DROP_INDEL_NOT_ON_REFERENCE))
        .when(untrusted_palindromic_mnp)
        .then(pl.lit(DROP_PALINDROMIC_MNP_UNTRUSTED))
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias(DROP_REASON_COL)
    )


def _apply_allele_actions(frame: pl.DataFrame, rules: Mapping[str, FlipRule]) -> pl.DataFrame:
    ea = pl.col(GWASLAB_EFFECT_ALLELE_COL)
    nea = pl.col(GWASLAB_NON_EFFECT_ALLELE_COL)
    complement = pl.col(ACTION_COL).is_in([ACTION_COMPLEMENT, ACTION_COMPLEMENT_SWAP])
    swap = pl.col(ACTION_COL).is_in([ACTION_SWAP, ACTION_COMPLEMENT_SWAP])
    complemented = frame.with_columns(
        pl.when(complement)
        .then(reverse_complement_expr(GWASLAB_EFFECT_ALLELE_COL))
        .otherwise(ea)
        .alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(complement)
        .then(reverse_complement_expr(GWASLAB_NON_EFFECT_ALLELE_COL))
        .otherwise(nea)
        .alias(GWASLAB_NON_EFFECT_ALLELE_COL),
    )
    swapped = complemented.with_columns(
        pl.when(swap).then(nea).otherwise(ea).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(swap).then(ea).otherwise(nea).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
    )
    return flip_statistics(swapped, mask=swap, rules=rules)
```

- [ ] **Step 10: Implement genome_reference_harmonization_task.py**

```python
"""
Genome-reference harmonization of GWAS summary statistics.

Orients every variant so that NEA is the plus-strand reference allele of a genome
FASTA, flips allele-dependent statistics to match, and resolves the variants the genome
alone cannot orient.

- **Trust.** A table whose checkable SNVs and indels are 100% reference-consistent is
  trusted: its palindromic variants and its indels with both alleles on the genome keep
  their source orientation.
- **Untrusted tables.** Palindromic SNV strands and ambiguous indels are resolved against
  a reference panel's allele frequencies, and are dropped whenever the evidence is not
  decisive.

This replaces gwaslab harmonization. See the design spec in
experiments/claude/design_specs/2026-09-14-genome-reference-harmonization-design.md.

Memory is bounded by one chromosome. Pass 1 reads four columns per chromosome to decide
trust. Pass 2 resolves one chromosome at a time into parquet parts, which are then
concatenated with a streaming sink.
"""

from collections.abc import Sequence
from pathlib import Path

import attrs
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    classify_alleles,
    prepare_alleles,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    DROPPED_COLUMNS,
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
    resolve_chromosome,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import (
    TrustCounts,
    count_trust_evidence,
    decide_trust,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

logger = structlog.get_logger()

HARMONIZED_FILENAME = "harmonized.parquet"
_KEPT_LABEL = "kept"  # log label for rows without a drop reason
_KEY_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]


def scan_sumstats_as_polars(asset: Asset, meta: Meta, pipe: DataProcessingPipe) -> pl.LazyFrame:
    """Scan the sumstats asset, apply the pipe, and insist the result is still polars-backed."""
    native = pipe.process(scan_dataframe_asset(asset, meta)).to_native()
    assert isinstance(native, pl.LazyFrame), (
        "genome-reference harmonization streams polars LazyFrames, but the pipe produced "
        f"{type(native).__name__}"
    )
    return native


def chromosomes_to_harmonize(
    sumstats: pl.LazyFrame, fasta: IndexedFasta, options: GenomeReferenceHarmonizationOptions
) -> list[int]:
    null_counts = (
        sumstats.select([pl.col(column).null_count() for column in _KEY_COLUMNS])
        .collect(engine="streaming")
        .row(0, named=True)
    )
    assert all(count == 0 for count in null_counts.values()), (
        f"null values in key columns: {null_counts}"
    )
    present = sorted(
        int(chrom)
        for chrom in sumstats.select(pl.col(GWASLAB_CHROM_COL).unique())
        .collect(engine="streaming")[GWASLAB_CHROM_COL]
        .to_list()
    )
    excluded = [chrom for chrom in present if chrom in options.excluded_chromosomes]
    if excluded:
        logger.info("dropping rows on excluded chromosomes", chromosomes=excluded)
    kept = [chrom for chrom in present if chrom not in options.excluded_chromosomes]
    missing = [chrom for chrom in kept if chrom not in fasta.entries]
    assert not missing, f"chromosomes {missing} are not in the FASTA {fasta.fasta_path}"
    return kept


def count_trust_evidence_genome_wide(
    sumstats: pl.LazyFrame,
    chromosomes: Sequence[int],
    fasta: IndexedFasta,
    options: GenomeReferenceHarmonizationOptions,
) -> TrustCounts:
    total = TrustCounts.zero()
    for chrom in chromosomes:
        keys = (
            sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
            .select(_KEY_COLUMNS)
            .collect(engine="streaming")
        )
        valid = prepare_alleles(keys).filter(valid_alleles_expr())
        classified = classify_alleles(
            valid, fasta=fasta, chrom=chrom, max_gather_bytes=options.max_gather_bytes
        )
        total = total + count_trust_evidence(classified)
    return total


@frozen
class ParquetPanelLoader:
    """Reads one chromosome's panel rows at requested positions from the panel parquet."""

    panel_path: Path
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
                PANEL_AF_COL,
            )
            .collect()
        )


def resolve_chromosome_rows(
    sumstats: pl.LazyFrame, context: ChromosomeContext, panel_path: Path
) -> pl.DataFrame:
    """All rows of one chromosome, resolved, with DROP_REASON_COL (used by experiments too)."""
    rows = sumstats.filter(pl.col(GWASLAB_CHROM_COL) == context.chrom).collect(
        engine="streaming"
    )
    return resolve_chromosome(
        rows, context, load_panel=ParquetPanelLoader(panel_path=panel_path, chrom=context.chrom)
    )


def _write_chromosome_part(
    sumstats: pl.LazyFrame, context: ChromosomeContext, panel_path: Path, parts_dir: Path
) -> Path:
    resolved = resolve_chromosome_rows(sumstats, context, panel_path)
    rows_by_reason = {
        (_KEPT_LABEL if reason is None else reason): count
        for reason, count in resolved.group_by(DROP_REASON_COL).len().rows()
    }
    logger.info(
        "genome-reference harmonization chromosome summary",
        chromosome=context.chrom,
        rows_by_drop_reason=rows_by_reason,
    )
    kept = (
        resolved.filter(pl.col(DROP_REASON_COL).is_null())
        .drop(DROP_REASON_COL, *[column for column in DROPPED_COLUMNS if column in resolved.columns])
        .sort(GWASLAB_POS_COL)
    )
    duplicated = kept.select(
        pl.struct(GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL)
        .is_duplicated()
        .any()
    ).item()
    assert not duplicated, (
        f"chromosome {context.chrom}: (CHR, POS, EA, NEA) is not unique after harmonization"
    )
    part_path = parts_dir / f"chr{context.chrom}.parquet"
    kept.write_parquet(part_path)
    return part_path


@frozen
class GenomeReferenceHarmonizationTask(Task):
    meta: FilteredGWASDataMeta
    sumstats_task: Task
    fasta_task: Task
    panel_task: Task
    options: GenomeReferenceHarmonizationOptions
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def deps(self) -> list[Task]:
        return [self.sumstats_task, self.fasta_task, self.panel_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats = scan_sumstats_as_polars(
            fetch(self.sumstats_task.asset_id), self.sumstats_task.meta, self.pipe
        )
        fasta_asset = fetch(self.fasta_task.asset_id)
        assert isinstance(fasta_asset, DirectoryAsset), (
            f"expected {self.fasta_task.asset_id} to be a DirectoryAsset"
        )
        panel_asset = fetch(self.panel_task.asset_id)
        assert isinstance(panel_asset, FileAsset), (
            f"expected {self.panel_task.asset_id} to be a FileAsset"
        )
        fasta = IndexedFasta.open(fasta_asset.path)
        rules = resolve_column_rules(
            columns=sumstats.collect_schema().names(), extra=self.options.extra_column_rules
        )
        chromosomes = chromosomes_to_harmonize(sumstats, fasta, self.options)
        counts = count_trust_evidence_genome_wide(sumstats, chromosomes, fasta, self.options)
        trusted = decide_trust(counts, self.options)
        logger.info(
            "genome-reference harmonization trust decision",
            trusted=trusted,
            counts=attrs.asdict(counts),
        )
        parts_dir = scratch_dir / "parts"
        parts_dir.mkdir()
        part_paths = [
            _write_chromosome_part(
                sumstats,
                ChromosomeContext(
                    chrom=chrom, fasta=fasta, trusted=trusted, rules=rules, options=self.options
                ),
                panel_path=panel_asset.path,
                parts_dir=parts_dir,
            )
            for chrom in chromosomes
        ]
        assert part_paths, "no chromosomes to harmonize"
        out_path = scratch_dir / HARMONIZED_FILENAME
        pl.concat([pl.scan_parquet(path) for path in part_paths]).sink_parquet(out_path)
        n_rows = pl.scan_parquet(out_path).select(pl.len()).collect().item()
        assert n_rows > 0, "no variants survived genome-reference harmonization"
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sumstats_task: Task,
        fasta_task: Task,
        panel_task: Task,
        options: GenomeReferenceHarmonizationOptions = GenomeReferenceHarmonizationOptions(),
        pipe: DataProcessingPipe = IdentityPipe(),
    ) -> "GenomeReferenceHarmonizationTask":
        source_meta = sumstats_task.meta
        assert isinstance(source_meta, FilteredGWASDataMeta), (
            f"expected a FilteredGWASDataMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir="processed",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            sumstats_task=sumstats_task,
            fasta_task=fasta_task,
            panel_task=panel_task,
            options=options,
            pipe=pipe,
        )
```

- [ ] **Step 11: Run the tests to verify they pass**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization -q`

Expected: all pass (Task 1, 2 and 4 tests).

If a test fails, fix the implementation, not the expectation. The fixture
comments document why each expected value is right.

- [ ] **Step 12: Run the full check**

Run:

```bash
pixi r invoke green > /tmp/claude-green.log 2>&1; echo "EXIT=$?"
grep -aE "passed|failed|error|All checks passed" /tmp/claude-green.log | tail
```

Expected: EXIT=0 and a pytest summary with no failures. Fix lint, format,
spellcheck (typos) and ty findings in the new files.

- [ ] **Step 13: Commit**

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization test_mecfs_bio/unit/build_system/task/genome_reference_harmonization
git commit -m "Add streaming GenomeReferenceHarmonizationTask core"
```

### Task 5: Untrusted palindromic SNV strand rules

**Files:**
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/palindromes.py
- Modify: mecfs_bio/build_system/task/genome_reference_harmonization/resolve_chromosome.py
- Modify (append tests): test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py

**Interfaces:**
- Consumes: resolve_chromosome internals from Task 4 (ACTION_COL, DROP_REASON_COL, ROW_INDEX_COL, IS_PALINDROMIC_SNV_COL, flip_statistics); PanelLoader.
- Produces:
  - decide_palindrome_strands(rows: pl.DataFrame, panel: pl.DataFrame, options) -> pl.Series, aligned with rows. rows holds POS, EA, NEA, EAF, oriented so that NEA is the reference base.
  - In resolve_chromosome.py: a panel is loaded once per untrusted chromosome, before allele actions, via _panel_positions(valid). Task 6 reuses it.

- [ ] **Step 1: Append the failing tests**

Add PanelRecord to the fixtures import list at the top of the test module, then
append:

```python
_PALINDROME_PANEL = [
    PanelRecord(pos=13, ref="T", alt="A", af=0.15),
    PanelRecord(pos=14, ref="G", alt="C", af=0.85),
    PanelRecord(pos=18, ref="T", alt="A", af=0.1),
    PanelRecord(pos=19, ref="T", alt="A", af=0.5),
    PanelRecord(pos=27, ref="T", alt="A", af=0.45),
]
_PALINDROMES = [
    Variant(pos=13, ea="A", nea="T", eaf=0.1, beta=0.3),  # same side of 0.5 as the panel -> kept
    Variant(pos=14, ea="G", nea="C", eaf=0.8, beta=0.3),  # EA is ref: swapped, then strand-flipped
    Variant(pos=18, ea="A", nea="T", eaf=0.9, beta=0.3),  # opposite side of 0.5 -> strand-flipped
    Variant(pos=19, ea="A", nea="T", eaf=0.5),  # sumstats MAF above 0.4 -> unresolved
    Variant(pos=24, ea="G", nea="C", eaf=0.1),  # no panel record -> unresolved
    Variant(pos=27, ea="A", nea="T", eaf=0.1),  # panel MAF above 0.4 -> unresolved
    Variant(pos=30, ea="T", nea="A", eaf=None),  # no EAF -> unresolved
]
_UNRESOLVED_PALINDROME_POSITIONS = [19, 24, 27, 30]


def test_untrusted_palindromes_are_resolved_by_panel_frequency(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, *_PALINDROMES]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants), panel=_PALINDROME_PANEL)
    assert positions(result) == [1, 2, 13, 14, 18]
    kept = row_at(result, 13)
    assert (kept[EA], kept[NEA], kept[BETA]) == ("A", "T", pytest.approx(0.3))
    swapped_then_flipped = row_at(result, 14)
    assert (swapped_then_flipped[EA], swapped_then_flipped[NEA]) == ("C", "G")
    assert swapped_then_flipped[BETA] == pytest.approx(0.3)
    assert swapped_then_flipped[EAF] == pytest.approx(0.8)
    flipped = row_at(result, 18)
    assert (flipped[EA], flipped[NEA]) == ("A", "T")
    assert (flipped[BETA], flipped[EAF]) == (pytest.approx(-0.3), pytest.approx(0.1))


def test_unresolved_palindromes_can_be_kept(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, *_PALINDROMES]
    options = attrs.evolve(TEST_OPTIONS, keep_unresolved_palindromes=True)
    result = run_harmonization(
        tmp_path / "run", sumstats_frame(variants), panel=_PALINDROME_PANEL, options=options
    )
    assert set(_UNRESOLVED_PALINDROME_POSITIONS) <= set(positions(result))


def test_trusted_palindromes_keep_source_strand(tmp_path: Path) -> None:
    opposite_side = Variant(pos=18, ea="A", nea="T", eaf=0.9, beta=0.3)
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, opposite_side]),
        panel=_PALINDROME_PANEL,
    )
    assert row_at(result, 18)[BETA] == pytest.approx(0.3)
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py -q -k palindrome`

Expected:
- test_untrusted_palindromes_are_resolved_by_panel_frequency FAILS: positions include 19, 24, 27, 30, and 18 is not flipped.
- The other two pass already.

- [ ] **Step 3: Implement palindromes.py**

```python
"""
Strand resolution of palindromic SNVs in untrusted tables, following gwaslab's rule.

Rows arrive oriented so that NEA is the reference base, so a panel record with
REF = NEA and ALT = EA describes the same variant. A variant is resolvable only when:

- its EAF is at most palindrome_maf_threshold, or at least 1 - palindrome_maf_threshold;
- the panel record exists;
- the panel MAF is at most panel_maf_threshold.

It is kept when EAF and panel AF lie on the same side of 0.5, and strand-flipped
(statistics flipped, alleles unchanged) otherwise.
"""

import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    PALINDROME_KEEP,
    PALINDROME_STRAND_FLIP,
    PALINDROME_UNRESOLVED,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

FREQUENCY_EPSILON = 1e-6
_PANEL_AF = "_panel_af"
_DECISION = "decision"
_KEYS = [GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]


def decide_palindrome_strands(
    rows: pl.DataFrame, panel: pl.DataFrame, options: GenomeReferenceHarmonizationOptions
) -> pl.Series:
    """One PalindromeDecision per row of rows (POS, EA, NEA, EAF), in row order."""
    records = panel.select(
        GWASLAB_POS_COL,
        pl.col(PANEL_REF_COL).alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.col(PANEL_ALT_COL).alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.col(PANEL_AF_COL).cast(pl.Float64).alias(_PANEL_AF),
    )
    joined = rows.join(records, on=_KEYS, how="left", maintain_order="left")
    assert joined.height == rows.height, "duplicate panel records for a palindromic SNV"
    eaf = pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    af = pl.col(_PANEL_AF)
    threshold = options.palindrome_maf_threshold
    eaf_informative = (eaf <= threshold + FREQUENCY_EPSILON) | (
        eaf >= 1 - threshold - FREQUENCY_EPSILON
    )
    panel_informative = (
        pl.min_horizontal(af, 1 - af) <= options.panel_maf_threshold + FREQUENCY_EPSILON
    )
    resolvable = (eaf_informative & panel_informative).fill_null(False)
    same_side = ((af < 0.5) & (eaf < 0.5)) | ((af > 0.5) & (eaf > 0.5))
    return joined.select(
        pl.when(~resolvable)
        .then(pl.lit(PALINDROME_UNRESOLVED))
        .when(same_side)
        .then(pl.lit(PALINDROME_KEEP))
        .otherwise(pl.lit(PALINDROME_STRAND_FLIP))
        .alias(_DECISION)
    )[_DECISION]
```

- [ ] **Step 4: Wire the rules into resolve_chromosome.py**

Add imports:

```python
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    CLASS_INDEL_BOTH,
    IS_PALINDROMIC_SNV_COL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    DROP_PALINDROME_UNRESOLVED,
    PALINDROME_STRAND_FLIP,
    PALINDROME_UNRESOLVED,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.palindromes import (
    decide_palindrome_strands,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_POS_COL,
)
```

Add the constant:

```python
PALINDROME_DECISION_COL = "_palindrome_decision"
```

In resolve_chromosome, replace

```python
    valid = _apply_allele_actions(valid, context.rules)
```

with

```python
    panel = None if context.trusted else load_panel(_panel_positions(valid))
    valid = _apply_allele_actions(valid, context.rules)
    if panel is not None:
        valid = _apply_palindrome_strand_rules(valid, panel, context)
```

Add the helpers:

```python
def _panel_positions(valid: pl.DataFrame) -> pl.Series:
    """Positions whose variants may need panel evidence: palindromic SNVs and ambiguous indels."""
    needs_panel = pl.col(IS_PALINDROMIC_SNV_COL) | (pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH)
    return valid.filter(needs_panel)[GWASLAB_POS_COL]


def _eaf_expr(frame: pl.DataFrame) -> pl.Expr:
    if GWASLAB_EFFECT_ALLELE_FREQ_COL in frame.columns:
        return pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    return pl.lit(None, dtype=pl.Float64).alias(GWASLAB_EFFECT_ALLELE_FREQ_COL)


def _apply_palindrome_strand_rules(
    valid: pl.DataFrame, panel: pl.DataFrame, context: ChromosomeContext
) -> pl.DataFrame:
    candidates = valid.filter(pl.col(IS_PALINDROMIC_SNV_COL) & pl.col(DROP_REASON_COL).is_null())
    decisions = candidates.select(ROW_INDEX_COL).with_columns(
        decide_palindrome_strands(
            candidates.select(
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
                _eaf_expr(candidates),
            ),
            panel,
            context.options,
        ).alias(PALINDROME_DECISION_COL)
    )
    decided = valid.join(decisions, on=ROW_INDEX_COL, how="left", maintain_order="left")
    decided = flip_statistics(
        decided, mask=pl.col(PALINDROME_DECISION_COL) == PALINDROME_STRAND_FLIP, rules=context.rules
    )
    if not context.options.keep_unresolved_palindromes:
        unresolved = (pl.col(PALINDROME_DECISION_COL) == PALINDROME_UNRESOLVED) & pl.col(
            DROP_REASON_COL
        ).is_null()
        decided = decided.with_columns(
            pl.when(unresolved)
            .then(pl.lit(DROP_PALINDROME_UNRESOLVED))
            .otherwise(pl.col(DROP_REASON_COL))
            .alias(DROP_REASON_COL)
        )
    return decided.drop(PALINDROME_DECISION_COL)
```

Merge the new imports into the existing import blocks rather than duplicating
module imports; ruff will flag duplicates.

- [ ] **Step 5: Run the module tests**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization -q`
Expected: all pass.

- [ ] **Step 6: Commit**

```bash
git add mecfs_bio/build_system/task/genome_reference_harmonization test_mecfs_bio/unit/build_system/task/genome_reference_harmonization
git commit -m "Resolve untrusted palindromic SNV strands by panel frequency"
```

### Task 6: Stringent ambiguous-indel rules and trust-decision tests

**Files:**
- Create: mecfs_bio/build_system/task/genome_reference_harmonization/ambiguous_indels.py
- Modify: mecfs_bio/build_system/task/genome_reference_harmonization/resolve_chromosome.py
- Modify (append tests): test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py

**Interfaces:**
- Consumes: the panel loaded in resolve_chromosome (Task 5); _eaf_expr; ACTION_COL; DROP_REASON_COL.
- Produces:
  - INDEL_ACTION_COL = "_indel_action"
  - INDEL_DROP_REASON_COL = "_indel_drop_reason"
  - decide_ambiguous_indels(rows: pl.DataFrame, panel: pl.DataFrame, options) -> pl.DataFrame. Output is aligned with rows (POS, EA, NEA, EAF) and has exactly those two columns: action ACTION_KEEP, ACTION_SWAP or null, and the drop reason or null.

- [ ] **Step 1: Append the failing tests**

```python
_AMBIGUOUS_DISTANCE = 0.1
_AMBIGUOUS_MARGIN = 0.2
_STRINGENT_OPTIONS = attrs.evolve(
    TEST_OPTIONS,
    indel_max_af_distance=_AMBIGUOUS_DISTANCE,
    indel_min_af_margin=_AMBIGUOUS_MARGIN,
)
# Every variant below is T/TT inside the chr1 T homopolymer at 31-40, so both alleles
# match the genome. "keep" reads the row as REF=TT, ALT=T; "flip" as REF=T, ALT=TT.
_AMBIGUOUS_INDELS = [
    Variant(pos=31, ea="T", nea="TT", eaf=None),  # no EAF -> dropped
    Variant(pos=32, ea="T", nea="TT", eaf=0.3),  # no panel record -> dropped
    Variant(pos=33, ea="T", nea="TT", eaf=0.3),  # keep record fits -> kept
    Variant(pos=34, ea="T", nea="TT", eaf=0.3, beta=0.2),  # flip record fits -> swapped
    Variant(pos=35, ea="T", nea="TT", eaf=0.3),  # only record misfits -> dropped
    Variant(pos=36, ea="T", nea="TT", eaf=0.3),  # both records, keep decisive -> kept
    Variant(pos=37, ea="T", nea="TT", eaf=0.5),  # both records, margin too small -> dropped
]
_AMBIGUOUS_PANEL = [
    PanelRecord(pos=33, ref="TT", alt="T", af=0.32),
    PanelRecord(pos=34, ref="T", alt="TT", af=0.68),
    PanelRecord(pos=35, ref="TT", alt="T", af=0.8),
    PanelRecord(pos=36, ref="TT", alt="T", af=0.31),
    PanelRecord(pos=36, ref="T", alt="TT", af=0.0),
    PanelRecord(pos=37, ref="TT", alt="T", af=0.5),
    PanelRecord(pos=37, ref="T", alt="TT", af=0.45),
]


def test_untrusted_ambiguous_indels_follow_the_stringent_rules(tmp_path: Path) -> None:
    variants = [CONSISTENT_SNV, INCONSISTENT_SNV, *_AMBIGUOUS_INDELS]
    result = run_harmonization(
        tmp_path / "run", sumstats_frame(variants), panel=_AMBIGUOUS_PANEL, options=_STRINGENT_OPTIONS
    )
    assert positions(result) == [1, 2, 33, 34, 36]
    swapped = row_at(result, 34)
    assert (swapped[EA], swapped[NEA]) == ("TT", "T")
    assert (swapped[EAF], swapped[BETA]) == (pytest.approx(0.7), pytest.approx(-0.2))


def test_ambiguous_indel_resolution_is_symmetric_in_orientation(tmp_path: Path) -> None:
    as_deletion_label = Variant(pos=36, ea="T", nea="TT", eaf=0.3, beta=0.2)
    as_insertion_label = Variant(pos=36, ea="TT", nea="T", eaf=0.7, beta=-0.2)
    results = [
        row_at(
            run_harmonization(
                tmp_path / name,
                sumstats_frame([CONSISTENT_SNV, INCONSISTENT_SNV, variant]),
                panel=_AMBIGUOUS_PANEL,
                options=_STRINGENT_OPTIONS,
            ),
            36,
        )
        for name, variant in [("deletion", as_deletion_label), ("insertion", as_insertion_label)]
    ]
    for row in results:
        assert (row[EA], row[NEA]) == ("T", "TT")
        assert (row[EAF], row[BETA]) == (pytest.approx(0.3), pytest.approx(0.2))


_AMBIGUOUS_WITHOUT_PANEL = Variant(pos=5, ea="T", nea="TG")


@pytest.mark.parametrize(
    ("extra_variants", "options", "expect_trusted"),
    [
        ([], TEST_OPTIONS, True),
        ([INCONSISTENT_SNV], TEST_OPTIONS, False),
        ([], attrs.evolve(TEST_OPTIONS, min_checkable_snvs=2), False),
        ([], attrs.evolve(TEST_OPTIONS, min_checkable_indels=2), False),
    ],
)
def test_trust_decision_controls_ambiguous_indels(
    tmp_path: Path,
    extra_variants: list[Variant],
    options: GenomeReferenceHarmonizationOptions,
    expect_trusted: bool,
) -> None:
    variants = [CONSISTENT_SNV, CONSISTENT_INDEL, _AMBIGUOUS_WITHOUT_PANEL, *extra_variants]
    result = run_harmonization(tmp_path / "run", sumstats_frame(variants), options=options)
    # Trusted tables keep the ambiguous indel; untrusted ones drop it (no panel record).
    assert (5 in positions(result)) == expect_trusted
```

Add this import at the top of the test module:

```python
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
```

- [ ] **Step 2: Run the new tests to verify they fail**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization/test_genome_reference_harmonization_task.py -q -k "ambiguous or trust_decision"`

Expected:
- the stringent-rules test FAILS: all T/TT rows are kept;
- the trust-decision cases with expect_trusted False FAIL.

- [ ] **Step 3: Implement ambiguous_indels.py**

```python
"""
Stringent resolution of ambiguous indels (both alleles on the genome) in untrusted tables.

Two readings of a row (EA, NEA, EAF) are possible.

- keep: the variant is REF=NEA, ALT=EA. Its panel record predicts EAF = AF.
- flip: the variant is REF=EA, ALT=NEA. Its panel record predicts EAF = 1 - AF, since EA
  is then the reference allele.

A reading is chosen only if its distance |EAF - predicted| is at most
indel_max_af_distance and, when the other reading also has a record, beats it by at
least indel_min_af_margin. Anything else is dropped: excluding correct variants is
preferred to keeping one whose orientation is wrong.
"""

import polars as pl

from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_KEEP,
    ACTION_SWAP,
    DROP_AMBIGUOUS_INDEL_AF_INDECISIVE,
    DROP_AMBIGUOUS_INDEL_AF_MISMATCH,
    DROP_AMBIGUOUS_INDEL_NO_EAF,
    DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

INDEL_ACTION_COL = "_indel_action"
INDEL_DROP_REASON_COL = "_indel_drop_reason"
_AF_KEEP = "_af_keep"
_AF_FLIP = "_af_flip"
_KEYS = [GWASLAB_POS_COL, GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL]


def _records(panel: pl.DataFrame, ref_as: str, alt_as: str, af_name: str) -> pl.DataFrame:
    return panel.select(
        GWASLAB_POS_COL,
        pl.col(PANEL_REF_COL).alias(ref_as),
        pl.col(PANEL_ALT_COL).alias(alt_as),
        pl.col(PANEL_AF_COL).cast(pl.Float64).alias(af_name),
    )


def decide_ambiguous_indels(
    rows: pl.DataFrame, panel: pl.DataFrame, options: GenomeReferenceHarmonizationOptions
) -> pl.DataFrame:
    """INDEL_ACTION_COL and INDEL_DROP_REASON_COL for each row of rows (POS, EA, NEA, EAF), in order."""
    keep_records = _records(
        panel, ref_as=GWASLAB_NON_EFFECT_ALLELE_COL, alt_as=GWASLAB_EFFECT_ALLELE_COL, af_name=_AF_KEEP
    )
    flip_records = _records(
        panel, ref_as=GWASLAB_EFFECT_ALLELE_COL, alt_as=GWASLAB_NON_EFFECT_ALLELE_COL, af_name=_AF_FLIP
    )
    joined = rows.join(keep_records, on=_KEYS, how="left", maintain_order="left").join(
        flip_records, on=_KEYS, how="left", maintain_order="left"
    )
    assert joined.height == rows.height, "duplicate panel records for an ambiguous indel"
    eaf = pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64)
    has_keep = pl.col(_AF_KEEP).is_not_null()
    has_flip = pl.col(_AF_FLIP).is_not_null()
    keep_distance = (eaf - pl.col(_AF_KEEP)).abs()
    flip_distance = (eaf - (1 - pl.col(_AF_FLIP))).abs()
    tolerance = options.indel_max_af_distance
    margin = options.indel_min_af_margin
    choose_keep = (
        has_keep
        & (keep_distance <= tolerance)
        & (~has_flip | (flip_distance - keep_distance >= margin))
    ).fill_null(False)
    choose_flip = (
        has_flip
        & (flip_distance <= tolerance)
        & (~has_keep | (keep_distance - flip_distance >= margin))
    ).fill_null(False)
    return joined.select(
        pl.when(choose_keep)
        .then(pl.lit(ACTION_KEEP))
        .when(choose_flip)
        .then(pl.lit(ACTION_SWAP))
        .otherwise(pl.lit(None, dtype=pl.String))
        .alias(INDEL_ACTION_COL),
        pl.when(eaf.is_null())
        .then(pl.lit(DROP_AMBIGUOUS_INDEL_NO_EAF))
        .when(~has_keep & ~has_flip)
        .then(pl.lit(DROP_AMBIGUOUS_INDEL_NOT_IN_PANEL))
        .when(choose_keep | choose_flip)
        .then(pl.lit(None, dtype=pl.String))
        .when(has_keep ^ has_flip)
        .then(pl.lit(DROP_AMBIGUOUS_INDEL_AF_MISMATCH))
        .otherwise(pl.lit(DROP_AMBIGUOUS_INDEL_AF_INDECISIVE))
        .alias(INDEL_DROP_REASON_COL),
    )
```

- [ ] **Step 4: Wire the rules into resolve_chromosome.py**

Add the import:

```python
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    INDEL_DROP_REASON_COL,
    decide_ambiguous_indels,
)
```

In resolve_chromosome, replace

```python
    valid = _apply_allele_actions(valid, context.rules)
```

with

```python
    if panel is not None:
        valid = _apply_ambiguous_indel_rules(valid, panel, context)
    valid = _apply_allele_actions(valid, context.rules)
```

Add the helper:

```python
def _apply_ambiguous_indel_rules(
    valid: pl.DataFrame, panel: pl.DataFrame, context: ChromosomeContext
) -> pl.DataFrame:
    ambiguous = valid.filter(
        (pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH) & pl.col(DROP_REASON_COL).is_null()
    )
    decisions = pl.concat(
        [
            ambiguous.select(ROW_INDEX_COL),
            decide_ambiguous_indels(
                ambiguous.select(
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                    _eaf_expr(ambiguous),
                ),
                panel,
                context.options,
            ),
        ],
        how="horizontal",
    )
    decided = valid.join(decisions, on=ROW_INDEX_COL, how="left", maintain_order="left")
    return decided.with_columns(
        pl.coalesce(pl.col(INDEL_ACTION_COL), pl.col(ACTION_COL)).alias(ACTION_COL),
        pl.coalesce(pl.col(DROP_REASON_COL), pl.col(INDEL_DROP_REASON_COL)).alias(DROP_REASON_COL),
    ).drop(INDEL_ACTION_COL, INDEL_DROP_REASON_COL)
```

- [ ] **Step 5: Run the module tests**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization -q`
Expected: all pass.

- [ ] **Step 6: Run green and commit**

```bash
pixi r invoke green > /tmp/claude-green.log 2>&1; echo "EXIT=$?"
grep -aE "passed|failed|error|All checks passed" /tmp/claude-green.log | tail
git add mecfs_bio/build_system/task/genome_reference_harmonization test_mecfs_bio/unit/build_system/task/genome_reference_harmonization
git commit -m "Apply stringent panel-frequency rules to untrusted ambiguous indels"
```

### Task 7: DropIndelsPipe

This replaces gwaslab filter_indels for callers of the annovar generator that
set filter_indels_in_harmonized.

**Files:**
- Create: mecfs_bio/build_system/task/pipes/drop_indels_pipe.py
- Create: test_mecfs_bio/unit/build_system/task/pipes/__init__.py (only if the directory does not exist)
- Create: test_mecfs_bio/unit/build_system/task/pipes/test_drop_indels_pipe.py

**Interfaces:**
- Produces: DropIndelsPipe(), a DataProcessingPipe that keeps rows whose EA and NEA are both one character.

- [ ] **Step 1: Write the failing test**

```python
import polars as pl

from mecfs_bio.build_system.task.pipes.drop_indels_pipe import DropIndelsPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)


def test_drop_indels_pipe_keeps_single_base_variants_only() -> None:
    frame = pl.DataFrame(
        {
            GWASLAB_POS_COL: [1, 2, 3, 4],
            GWASLAB_EFFECT_ALLELE_COL: ["A", "AT", "G", "GC"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C", "A", "GT", "TA"],
        }
    ).with_columns(pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.Categorical))
    result = DropIndelsPipe().process_eager_polars(frame)
    assert result[GWASLAB_POS_COL].to_list() == [1]
```

- [ ] **Step 2: Run it to verify it fails**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/pipes/test_drop_indels_pipe.py -q`
Expected: FAIL with ModuleNotFoundError.

- [ ] **Step 3: Implement the pipe**

```python
"""Drop variants whose effect or non-effect allele is longer than one base."""

import narwhals
from attrs import frozen

from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
)


@frozen
class DropIndelsPipe(DataProcessingPipe):
    """
    Keep only variants whose EA and NEA are both a single base.

    Useful ahead of genome-reference harmonization for datasets with very long
    structural-variant alleles that downstream SNP-based analyses do not use.
    """

    def process(self, x: narwhals.LazyFrame) -> narwhals.LazyFrame:
        return x.filter(
            (narwhals.col(GWASLAB_EFFECT_ALLELE_COL).cast(narwhals.String).str.len_chars() == 1)
            & (
                narwhals.col(GWASLAB_NON_EFFECT_ALLELE_COL)
                .cast(narwhals.String)
                .str.len_chars()
                == 1
            )
        )
```

- [ ] **Step 4: Run the test to verify it passes, then commit**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/pipes/test_drop_indels_pipe.py -q`
Expected: 1 passed.

```bash
git add mecfs_bio/build_system/task/pipes/drop_indels_pipe.py test_mecfs_bio/unit/build_system/task/pipes
git commit -m "Add DropIndelsPipe"
```

### Task 8: Validation V1 and V2, row-by-row parity with gwaslab harmonization

Run this before any call site is switched (Task 11): it needs the existing
gwaslab-harmonized assets and task objects.

**Files:**
- Create: experiments/claude/genome_reference_harmonization/compare_with_gwaslab.py
- Create: experiments/claude/genome_reference_harmonization/validation_report.md (written from the logs)

**Interfaces:**
- Consumes:
  - GenomeReferenceHarmonizationTask, scan_sumstats_as_polars, chromosomes_to_harmonize, count_trust_evidence_genome_wide, resolve_chromosome_rows, decide_trust, resolve_column_rules, ChromosomeContext, DROP_REASON_COL, reverse_complement_expr (Tasks 4-6)
  - the Task 3 assets
  - DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37, DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.dump_parquet_task
  - LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS, LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET

- [ ] **Step 1: Write the comparison script**

```python
"""
V1/V2: compare genome-reference harmonization with gwaslab harmonization, row by row.

For DecodeME (liftover to 37) and Liu et al. 2023 IBD:
- build the pre-harmonization table and run GenomeReferenceHarmonizationTask on it;
- recompute per-chromosome drop reasons with the library;
- join the output with the cached gwaslab-harmonized table on SNPID;
- bucket every row and break the buckets down by variant type.

The cached gwaslab tables may predate the installed gwaslab version, because the build
cache does not key on code.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.compare_with_gwaslab \
    2>&1 | tee experiments/claude/genome_reference_harmonization/compare_with_gwaslab.log
"""

from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_dump_to_parquet import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
)
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats import (
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_liftover_to_37 import (
    DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    reverse_complement_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import IndexedFasta
from mecfs_bio.build_system.task.genome_reference_harmonization.flip import (
    resolve_column_rules,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    resolve_chromosome_rows,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.resolve_chromosome import (
    DROP_REASON_COL,
    ChromosomeContext,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import decide_trust
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SNPID_COL,
    GWASLAB_STATUS_COL,
)
from mecfs_bio.constants.regenie_constants import (
    REGENIE_A1FREQ_CASES_COL,
    REGENIE_A1FREQ_CONTROLS_COL,
)

EA = GWASLAB_EFFECT_ALLELE_COL
NEA = GWASLAB_NON_EFFECT_ALLELE_COL
BETA = GWASLAB_BETA_COL
SUFFIX = "_gwaslab"
STAT_TOLERANCE = 1e-6
BUCKET_COL = "bucket"
VARIANT_TYPE_COL = "variant_type"
INDEL_INFERENCE_FLIP_COL = "gwaslab_indel_inference_flip"
TYPE_INDEL = "indel"
TYPE_PALINDROMIC_SNV = "palindromic_snv"
TYPE_SNV = "snv"
TYPE_MNP = "mnp"
BUCKET_NEW_ONLY = "new_only"
BUCKET_GWASLAB_ONLY = "gwaslab_only"
BUCKET_IDENTICAL = "identical"
BUCKET_SAME_ORIENTATION_BETA_DIFFERS = "same_orientation_beta_differs"
BUCKET_OPPOSITE_ORIENTATION = "opposite_orientation"
BUCKET_OTHER = "other"
# gwaslab STATUS digit 7 value 4: stats flipped by indel inference (the bug signature)
GWASLAB_STATUS_DIGIT7_INDEL_FLIPPED = "4"
OPTIONS = GenomeReferenceHarmonizationOptions()


@frozen
class ComparisonCase:
    label: str
    pre_harmonization: Task
    gwaslab_harmonized: Task


@frozen
class DropReport:
    trusted: bool
    dropped: pl.DataFrame


CASES = [
    ComparisonCase(
        label="decode_me_gwas_1",
        pre_harmonization=GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=DECODE_ME_GWAS_1_SUMSTATS_LIFTOVER_TO_37,
            asset_id="experiment_decode_me_gwas_1_pre_harmonization_table",
            sub_dir="processed",
        ),
        gwaslab_harmonized=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.dump_parquet_task,
    ),
    ComparisonCase(
        label="liu_et_al_2023_ibd",
        pre_harmonization=GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
            asset_id="experiment_liu_et_al_2023_ibd_pre_harmonization_table",
            sub_dir="processed",
        ),
        gwaslab_harmonized=LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
    ),
]


def _file(asset: Asset) -> Path:
    assert isinstance(asset, FileAsset)
    return asset.path


def _variant_type(ea: str, nea: str) -> pl.Expr:
    ea_length = pl.col(ea).str.len_bytes()
    nea_length = pl.col(nea).str.len_bytes()
    palindromic = (ea_length == 1) & (nea_length == 1) & (reverse_complement_expr(nea) == pl.col(ea))
    return (
        pl.when(ea_length != nea_length)
        .then(pl.lit(TYPE_INDEL))
        .when(palindromic)
        .then(pl.lit(TYPE_PALINDROMIC_SNV))
        .when(ea_length == 1)
        .then(pl.lit(TYPE_SNV))
        .otherwise(pl.lit(TYPE_MNP))
    )


def drop_report(case: ComparisonCase, assets: dict[str, Asset]) -> DropReport:
    pre_asset = assets[case.pre_harmonization.asset_id]
    fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    panel_path = _file(assets[THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES.asset_id])
    sumstats = scan_sumstats_as_polars(pre_asset, case.pre_harmonization.meta, IdentityPipe())
    fasta = IndexedFasta.open(fasta_asset.path)
    rules = resolve_column_rules(columns=sumstats.collect_schema().names(), extra=())
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
    counts = count_trust_evidence_genome_wide(sumstats, chromosomes, fasta, OPTIONS)
    trusted = decide_trust(counts, OPTIONS)
    print(f"[{case.label}] trust counts {counts} -> trusted={trusted}")
    parts = []
    for chrom in chromosomes:
        context = ChromosomeContext(chrom=chrom, fasta=fasta, trusted=trusted, rules=rules, options=OPTIONS)
        resolved = resolve_chromosome_rows(sumstats, context, panel_path)
        parts.append(
            resolved.filter(pl.col(DROP_REASON_COL).is_not_null()).select(
                GWASLAB_SNPID_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL, EA, NEA, DROP_REASON_COL
            )
        )
    return DropReport(trusted=trusted, dropped=pl.concat(parts))


def compare(case: ComparisonCase, new: pl.DataFrame, old: pl.DataFrame, drops: DropReport) -> None:
    for name, frame in [("new", new), ("gwaslab", old)]:
        assert frame[GWASLAB_SNPID_COL].is_unique().all(), f"{name} SNPIDs are not unique"
    old = old.with_columns(pl.col(EA).cast(pl.String), pl.col(NEA).cast(pl.String))
    joined = new.join(old, on=GWASLAB_SNPID_COL, how="full", suffix=SUFFIX, coalesce=True)
    same_orientation = (pl.col(EA) == pl.col(EA + SUFFIX)) & (pl.col(NEA) == pl.col(NEA + SUFFIX))
    opposite_orientation = (pl.col(EA) == pl.col(NEA + SUFFIX)) & (pl.col(NEA) == pl.col(EA + SUFFIX))
    bucket = (
        pl.when(pl.col(EA + SUFFIX).is_null())
        .then(pl.lit(BUCKET_NEW_ONLY))
        .when(pl.col(EA).is_null())
        .then(pl.lit(BUCKET_GWASLAB_ONLY))
        .when(same_orientation & ((pl.col(BETA) - pl.col(BETA + SUFFIX)).abs() < STAT_TOLERANCE))
        .then(pl.lit(BUCKET_IDENTICAL))
        .when(same_orientation)
        .then(pl.lit(BUCKET_SAME_ORIENTATION_BETA_DIFFERS))
        .when(opposite_orientation & ((pl.col(BETA) + pl.col(BETA + SUFFIX)).abs() < STAT_TOLERANCE))
        .then(pl.lit(BUCKET_OPPOSITE_ORIENTATION))
        .otherwise(pl.lit(BUCKET_OTHER))
    )
    joined = joined.with_columns(
        bucket.alias(BUCKET_COL),
        pl.coalesce(_variant_type(EA, NEA), _variant_type(EA + SUFFIX, NEA + SUFFIX)).alias(VARIANT_TYPE_COL),
        (pl.col(GWASLAB_STATUS_COL + SUFFIX).cast(pl.String).str.slice(6, 1) == GWASLAB_STATUS_DIGIT7_INDEL_FLIPPED).alias(
            INDEL_INFERENCE_FLIP_COL
        ),
    )
    with pl.Config(tbl_rows=60, tbl_cols=-1, tbl_width_chars=250):
        print(f"\n=== {case.label}: buckets by variant type ===")
        print(joined.group_by(BUCKET_COL, VARIANT_TYPE_COL).len().sort(BUCKET_COL, VARIANT_TYPE_COL))
        print(f"\n=== {case.label}: opposite-orientation rows by gwaslab indel-inference flip ===")
        print(
            joined.filter(pl.col(BUCKET_COL) == BUCKET_OPPOSITE_ORIENTATION)
            .group_by(VARIANT_TYPE_COL, INDEL_INFERENCE_FLIP_COL)
            .len()
        )
        gwaslab_only = joined.filter(pl.col(BUCKET_COL) == BUCKET_GWASLAB_ONLY).join(
            drops.dropped.select(GWASLAB_SNPID_COL, DROP_REASON_COL), on=GWASLAB_SNPID_COL, how="left"
        )
        print(f"\n=== {case.label}: rows only gwaslab kept, by our drop reason ===")
        print(gwaslab_only.group_by(DROP_REASON_COL, VARIANT_TYPE_COL).len().sort("len", descending=True))
        print(f"\n=== {case.label}: rows only we kept ===")
        print(joined.filter(pl.col(BUCKET_COL) == BUCKET_NEW_ONLY).group_by(VARIANT_TYPE_COL).len())
        for name in [BUCKET_SAME_ORIENTATION_BETA_DIFFERS, BUCKET_OTHER, BUCKET_NEW_ONLY, BUCKET_GWASLAB_ONLY]:
            print(f"\n--- {case.label}: examples of {name} ---")
            print(joined.filter(pl.col(BUCKET_COL) == name).head(8))
        for column in [REGENIE_A1FREQ_CASES_COL, REGENIE_A1FREQ_CONTROLS_COL]:
            if column in new.columns:
                differs = joined.filter(
                    (pl.col(BUCKET_COL) == BUCKET_IDENTICAL)
                    & ((pl.col(column) - pl.col(column + SUFFIX)).abs() > STAT_TOLERANCE)
                )
                print(f"\n{case.label}: identical rows whose {column} differs (gwaslab never flips it): {differs.height:,}")


def main() -> None:
    for case in CASES:
        new_task = GenomeReferenceHarmonizationTask.create(
            asset_id=f"experiment_{case.label}_genome_reference_harmonized",
            sumstats_task=case.pre_harmonization,
            fasta_task=UCSC_HG19_INDEXED_FASTA,
            panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
            options=OPTIONS,
        )
        assets = dict(
            DEFAULT_RUNNER.run(
                [
                    new_task,
                    case.pre_harmonization,
                    case.gwaslab_harmonized,
                    UCSC_HG19_INDEXED_FASTA,
                    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
                ]
            )
        )
        new = pl.read_parquet(_file(assets[new_task.asset_id]))
        old = pl.read_parquet(_file(assets[case.gwaslab_harmonized.asset_id]))
        compare(case, new, old, drop_report(case, assets))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the comparison**

Run it in the background; DecodeME takes minutes.

Run: `pixi r python -m experiments.claude.genome_reference_harmonization.compare_with_gwaslab 2>&1 | tee experiments/claude/genome_reference_harmonization/compare_with_gwaslab.log`

Expected:
- DecodeME trusted=False (spec: 99.83% SNV consistency).
- Buckets beyond "identical" are dominated by:
  - indel rows (opposite orientation where the gwaslab indel-inference flip is True; gwaslab_only with ambiguous_indel_* or indel_not_on_reference reasons);
  - palindromic SNVs;
  - A1FREQ differences.

If SNPIDs are not unique in either table, change the join key to (CHR, POS,
SNPID) and re-run.

- [ ] **Step 3: Write validation_report.md**

Record, per dataset:
- the trust counts and decision;
- the bucket and variant-type table;
- an explanation of every non-identical bucket, each tied to a spec rule.

The expected explanations are:
- the gwaslab indel-inference bug;
- stringent drops;
- reverse-strand indels;
- A1FREQ flips;
- invalid alleles;
- unknown contigs.

**Stop condition:** if any bucket cannot be explained by a spec rule (for
example non-indel opposite-orientation rows, or same_orientation_beta_differs),
stop and report it to the user before continuing. It is either a missed gwaslab
behaviour or a bug.

- [ ] **Step 4: Commit**

```bash
git add experiments/claude/genome_reference_harmonization
git commit -m "Validate genome-reference harmonization against gwaslab (V1, V2)"
```

### Task 9: Validation V3, tuning the stringent indel rules on a truth set

**Files:**
- Create: experiments/claude/genome_reference_harmonization/tune_ambiguous_indel_rules.py
- Modify: mecfs_bio/build_system/task/genome_reference_harmonization/options.py (defaults)
- Modify: experiments/claude/design_specs/2026-09-14-genome-reference-harmonization-design.md (record the chosen defaults)
- Modify: experiments/claude/genome_reference_harmonization/validation_report.md

**Interfaces:**
- Consumes:
  - classify_alleles, prepare_alleles, valid_alleles_expr, decide_ambiguous_indels, ParquetPanelLoader, chromosomes_to_harmonize, count_trust_evidence_genome_wide, decide_trust, scan_sumstats_as_polars (Tasks 4-6)
  - UCSC_HG38_INDEXED_FASTA, THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES (Task 3)
  - DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING

- [ ] **Step 1: Write the tuning script**

```python
"""
V3: tune indel_max_af_distance and indel_min_af_margin on a truth set.

DecodeME build 38 is 100% reference-consistent against hg38, so each ambiguous indel's
source orientation is the truth. Every ambiguous indel is resolved with the stringent
rules as if the table were untrusted, over a grid of options. Choosing "swap" is a wrong
choice. The chosen defaults are the grid cell with zero wrong choices that keeps the most
indels; ties go to the larger margin, then the smaller distance.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules \
    2>&1 | tee experiments/claude/genome_reference_harmonization/tune_ambiguous_indel_rules.log
"""

from pathlib import Path

import attrs
import polars as pl
from attrs import frozen

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_gwas_1_sumstats_minimal_processing import (
    DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_BOTH,
    classify_alleles,
    prepare_alleles,
    valid_alleles_expr,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    INDEL_DROP_REASON_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import IndexedFasta
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    ParquetPanelLoader,
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_KEEP,
    ACTION_SWAP,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import decide_trust
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

GRID_DISTANCES = (0.02, 0.05, 0.1, 0.15, 0.2)
GRID_MARGINS = (0.02, 0.05, 0.1, 0.2, 0.3)
BASE_OPTIONS = GenomeReferenceHarmonizationOptions()
ROW_COLUMNS = [
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
]

BUILD_38_TABLE = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=DECODE_ME_GWAS_1_SUMSTATS_MINIMAL_FILTERING,
    asset_id="experiment_decode_me_gwas_1_build_38_table",
    sub_dir="processed",
)


@frozen
class ChromosomeAmbiguousIndels:
    rows: pl.DataFrame
    panel: pl.DataFrame


def _ambiguous_indels(
    sumstats: pl.LazyFrame, chrom: int, fasta: IndexedFasta, panel_path: Path
) -> ChromosomeAmbiguousIndels:
    rows = (
        sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
        .select(ROW_COLUMNS)
        .collect(engine="streaming")
    )
    classified = classify_alleles(
        prepare_alleles(rows).filter(valid_alleles_expr()),
        fasta=fasta,
        chrom=chrom,
        max_gather_bytes=BASE_OPTIONS.max_gather_bytes,
    )
    ambiguous = classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH).select(ROW_COLUMNS)
    panel = ParquetPanelLoader(panel_path=panel_path, chrom=chrom)(ambiguous[GWASLAB_POS_COL])
    return ChromosomeAmbiguousIndels(rows=ambiguous, panel=panel)


def main() -> None:
    assets = DEFAULT_RUNNER.run(
        [BUILD_38_TABLE, UCSC_HG38_INDEXED_FASTA, THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES]
    )
    fasta_asset = assets[UCSC_HG38_INDEXED_FASTA.asset_id]
    panel_asset = assets[THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset) and isinstance(panel_asset, FileAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    sumstats = scan_sumstats_as_polars(
        assets[BUILD_38_TABLE.asset_id], BUILD_38_TABLE.meta, IdentityPipe()
    )
    chromosomes = chromosomes_to_harmonize(sumstats, fasta, BASE_OPTIONS)
    counts = count_trust_evidence_genome_wide(sumstats, chromosomes, fasta, BASE_OPTIONS)
    print(f"build-38 trust counts: {counts}")
    assert decide_trust(counts, BASE_OPTIONS), "build 38 is not trusted, so it is not a truth set"
    per_chromosome = [
        _ambiguous_indels(sumstats, chrom, fasta, panel_asset.path) for chrom in chromosomes
    ]
    n_ambiguous = sum(item.rows.height for item in per_chromosome)
    print(f"ambiguous indels: {n_ambiguous:,}")
    results = []
    for distance in GRID_DISTANCES:
        for margin in GRID_MARGINS:
            options = attrs.evolve(
                BASE_OPTIONS, indel_max_af_distance=distance, indel_min_af_margin=margin
            )
            decisions = pl.concat(
                [decide_ambiguous_indels(item.rows, item.panel, options) for item in per_chromosome]
            )
            reasons = dict(decisions.group_by(INDEL_DROP_REASON_COL).len().rows())
            results.append(
                {
                    "distance": distance,
                    "margin": margin,
                    "kept": int((decisions[INDEL_ACTION_COL] == ACTION_KEEP).sum()),
                    "wrong": int((decisions[INDEL_ACTION_COL] == ACTION_SWAP).sum()),
                    **{f"drop_{reason}": count for reason, count in reasons.items() if reason is not None},
                }
            )
    table = pl.DataFrame(results).fill_null(0)
    with pl.Config(tbl_rows=40, tbl_cols=-1, tbl_width_chars=250):
        print(table)
    candidates = table.filter(pl.col("wrong") == 0).sort(
        ["kept", "margin", "distance"], descending=[True, True, False]
    )
    assert candidates.height > 0, "no grid cell has zero wrong choices; report to the user"
    print("chosen defaults:", candidates.row(0, named=True))


if __name__ == "__main__":
    main()
```

- [ ] **Step 2: Run the tuning**

Run: `pixi r python -m experiments.claude.genome_reference_harmonization.tune_ambiguous_indel_rules 2>&1 | tee experiments/claude/genome_reference_harmonization/tune_ambiguous_indel_rules.log`

Expected:
- a trust line with zero inconsistent counts;
- the grid table;
- a "chosen defaults" line.

**Stop conditions** (stop and report to the user):
- the build-38 trust assertion fails;
- no cell has zero wrong choices.

- [ ] **Step 3: Set the defaults**

In options.py, set indel_max_af_distance and indel_min_af_margin to the chosen
cell's distance and margin.

Replace "Stringent indel options: indel_max_af_distance 0.1 and
indel_min_af_margin 0.2 until Task 9" in this plan's Global Constraints and the
"placeholders are 0.1 and 0.2" sentence in the design spec with the chosen
values and a pointer to the log.

Add the grid table and choice to validation_report.md.

- [ ] **Step 4: Run the module tests (they use explicit options) and commit**

Run: `pixi r python -m pytest test_mecfs_bio/unit/build_system/task/genome_reference_harmonization -q`
Expected: all pass.

```bash
git add experiments/claude/genome_reference_harmonization mecfs_bio/build_system/task/genome_reference_harmonization/options.py experiments/claude/design_specs
git commit -m "Tune stringent ambiguous-indel defaults on build-38 DecodeME (V3)"
```

### Task 10: Validation V4, peak memory

Run this before switching call sites; the gwaslab scenario uses the existing
gwaslab harmonization task object.

**Files:**
- Create: experiments/claude/genome_reference_harmonization/measure_harmonization_memory.py (child: one scenario, prints JSON)
- Create: experiments/claude/genome_reference_harmonization/memory_benchmark.py (driver)
- Modify: experiments/claude/genome_reference_harmonization/validation_report.md

**Interfaces:**
- Consumes: GenomeReferenceHarmonizationTask; the Task 8 pre-harmonization table (experiment_decode_me_gwas_1_pre_harmonization_table); DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task; FilterRowsByValue; execute_command.

- [ ] **Step 1: Write the child script**

```python
"""
Run one harmonization scenario in this process and print peak anonymous RSS as JSON.

RssAnon excludes file-backed pages, so the memory-mapped FASTA does not count against
the task. Scenarios:
- genome_reference_all: GenomeReferenceHarmonizationTask on DecodeME build 37
- genome_reference_chr1_2: the same, restricted to chromosomes 1 and 2 by a pipe
- gwaslab: the existing gwaslab harmonization task on the same input

Run (normally via memory_benchmark.py):
  pixi r python -m experiments.claude.genome_reference_harmonization.measure_harmonization_memory genome_reference_all
"""

import json
import sys
import tempfile
import threading
import time
from pathlib import Path

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.base_meta import DirMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL
from experiments.claude.genome_reference_harmonization.compare_with_gwaslab import CASES

SAMPLE_SECONDS = 0.2
_DECODE_ME_PRE_TABLE = CASES[0].pre_harmonization


def _rss_anon_kib() -> int:
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("RssAnon:"):
            return int(line.split()[1])
    raise RuntimeError("RssAnon not reported by /proc/self/status")


class _PeakSampler(threading.Thread):
    def __init__(self) -> None:
        super().__init__(daemon=True)
        self.peak_kib = 0
        self.stopped = threading.Event()

    def run(self) -> None:
        while not self.stopped.is_set():
            self.peak_kib = max(self.peak_kib, _rss_anon_kib())
            time.sleep(SAMPLE_SECONDS)


def _fetch_from_store(task: Task):
    by_id = {dep.asset_id: dep for dep in task.deps}

    def fetch(asset_id: AssetId) -> Asset:
        meta = by_id[asset_id].meta
        path = DEFAULT_RUNNER.meta_to_path(meta)
        return DirectoryAsset(path) if isinstance(meta, DirMeta) else FileAsset(path)

    return fetch


def _scenario_task(name: str) -> Task:
    if name == "gwaslab":
        return DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.harmonize_task
    pipe = (
        FilterRowsByValue(target_column=GWASLAB_CHROM_COL, valid_values=[1, 2])
        if name == "genome_reference_chr1_2"
        else IdentityPipe()
    )
    assert name in {"genome_reference_all", "genome_reference_chr1_2"}, f"unknown scenario {name}"
    return GenomeReferenceHarmonizationTask.create(
        asset_id=f"experiment_memory_{name}",
        sumstats_task=_DECODE_ME_PRE_TABLE,
        fasta_task=UCSC_HG19_INDEXED_FASTA,
        panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
        pipe=pipe,
    )


def main() -> None:
    name = sys.argv[1]
    task = _scenario_task(name)
    DEFAULT_RUNNER.run(list(task.deps))  # materialize inputs before sampling
    sampler = _PeakSampler()
    baseline_kib = _rss_anon_kib()
    sampler.start()
    started = time.monotonic()
    with tempfile.TemporaryDirectory() as scratch:
        task.execute(scratch_dir=Path(scratch), fetch=_fetch_from_store(task), wf=make_wf())
    elapsed = time.monotonic() - started
    sampler.stopped.set()
    sampler.join()
    print(
        json.dumps(
            {
                "scenario": name,
                "baseline_rss_anon_gib": baseline_kib / 2**20,
                "peak_rss_anon_gib": sampler.peak_kib / 2**20,
                "seconds": elapsed,
            }
        )
    )


if __name__ == "__main__":
    main()
```

DEFAULT_RUNNER.meta_to_path is the same accessor harmonizer_flip_audit.py
uses. If the runner exposes it under another name, use that; do not add a new
path helper.

- [ ] **Step 2: Write the driver**

```python
"""
V4: compare peak anonymous memory of genome-reference harmonization with gwaslab harmonization.

Each scenario runs in its own process. Acceptance: the all-chromosome peak is at most
1.5 times the chromosome-1-and-2 peak, i.e. peak memory follows the largest chromosome
rather than total rows.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.memory_benchmark \
    2>&1 | tee experiments/claude/genome_reference_harmonization/memory_benchmark.log
"""

import json

from mecfs_bio.util.subproc.run_command import execute_command

SCENARIOS = ["genome_reference_chr1_2", "genome_reference_all", "gwaslab"]
ACCEPTABLE_RATIO = 1.5


def main() -> None:
    peaks: dict[str, float] = {}
    for scenario in SCENARIOS:
        output = execute_command(
            [
                "python",
                "-m",
                "experiments.claude.genome_reference_harmonization.measure_harmonization_memory",
                scenario,
            ]
        )
        result = json.loads(output.strip().splitlines()[-1])
        print(result)
        peaks[scenario] = result["peak_rss_anon_gib"]
    ratio = peaks["genome_reference_all"] / peaks["genome_reference_chr1_2"]
    print(f"all / chr1-2 peak ratio: {ratio:.2f} (acceptable <= {ACCEPTABLE_RATIO})")
    print(f"gwaslab / genome-reference peak ratio: {peaks['gwaslab'] / peaks['genome_reference_all']:.2f}")
    assert ratio <= ACCEPTABLE_RATIO, "peak memory grows with total rows; investigate before switching"


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run the benchmark**

Run it in the background; the gwaslab scenario is slow.

Run: `pixi r python -m experiments.claude.genome_reference_harmonization.memory_benchmark 2>&1 | tee experiments/claude/genome_reference_harmonization/memory_benchmark.log`

Expected: three JSON lines, a ratio line, and no assertion error.

**Stop condition:** if the ratio assertion fails, stop and report to the user
with the log.

- [ ] **Step 4: Record and commit**

Add the three peaks, run times and the ratio to validation_report.md.

```bash
git add experiments/claude/genome_reference_harmonization
git commit -m "Measure genome-reference harmonization peak memory (V4)"
```

### Task 11: Switch call sites to genome-reference harmonization

Prerequisite: Tasks 8-10 finished with no stop condition triggered.

After this task, compare_with_gwaslab.py and measure_harmonization_memory.py no
longer import: they reference the removed gwaslab-harmonization task objects.
They are historical records of V1, V2 and V4, and their logs are the evidence;
do not update them. Add one line to each module docstring saying they ran
against the pre-switch code at the commit made in Task 10.

**Files:**
- Modify: mecfs_bio/asset_generator/annovar_37_basic_rsid_assignment.py
- Modify: mecfs_bio/assets/gwas/inflammatory_bowel_disease/liu_et_al_2023/processed_gwas_data/liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized.py
- Delete: mecfs_bio/assets/gwas/inflammatory_bowel_disease/liu_et_al_2023/processed_gwas_data/liu_et_al_2023_eur_37_harmonized_dump_to_parquet.py
- Modify: .../liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar.py and .../liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar_with_dups.py
- Modify: test_mecfs_bio/system/test_harmonize_drop_ambiguous.py
- Modify: .github/workflows/system_test_harmonize.yml

**Interfaces:**
- Consumes: GenomeReferenceHarmonizationTask, GenomeReferenceHarmonizationOptions, DropIndelsPipe, the Task 3 hg19 assets, IndexedFasta, reference_matches.
- Produces: RSIDAssignmentTaskGroup(pre_harmonization_table_task, harmonize_task, join_task). The dump_parquet_task field is replaced; grep confirmed no caller reads the group's harmonize_task or dump_parquet_task.

- [ ] **Step 1: Rewrite annovar_37_basic_rsid_assignment.py**

```python
"""
Asset generator for assigning rsIDs to genome build-37 GWAS datasets.
"""

import narwhals
from attrs import frozen

from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename_unique import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_DIRECT_DOWNLOAD,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.drop_indels_pipe import DropIndelsPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe


@frozen
class RSIDAssignmentTaskGroup:
    """
    Collection of tasks used to assign rsIDs by joining with an existing dataframe of SNPs
    """

    pre_harmonization_table_task: Task
    harmonize_task: Task
    join_task: Task


def annovar_37_basic_rsid_assignment(
    sumstats_task: Task,
    base_name: str,
    use_gwaslab_rsids_convention: bool = False,
    drop_palindromic_ambiguous: bool = True,
    filter_indels_in_harmonized: bool = False,
) -> RSIDAssignmentTaskGroup:
    """
    Asset generator that creates a chain of tasks to assign rsIDs to existing build 37
    sumstats datasets using the annovar dbSNP reference data.

    The gwaslab Sumstats object is dumped to a table and oriented by genome-reference
    harmonization against the UCSC hg19 FASTA and the 1000 Genomes EUR panel.

    Set drop_palindromic_ambiguous to False to keep palindromic SNVs whose strand cannot be
    resolved. Ambiguous indels are never kept on that basis.

    Set filter_indels_in_harmonized to drop indels before harmonization. This suits datasets
    with very long structural-variant alleles that the downstream SNP-based analyses (LDSC
    genetic correlation, MAGMA) do not use.
    """
    pre_harmonization_table_task = GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=sumstats_task,
        asset_id=base_name + "_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
        pipe=DropIndelsPipe() if filter_indels_in_harmonized else IdentityPipe(),
    )
    harmonize_task = GenomeReferenceHarmonizationTask.create(
        asset_id=base_name + "_genome_reference_harmonized",
        sumstats_task=pre_harmonization_table_task,
        fasta_task=UCSC_HG19_INDEXED_FASTA,
        panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
        options=GenomeReferenceHarmonizationOptions(
            keep_unresolved_palindromes=not drop_palindromic_ambiguous
        ),
    )
    out_pipe: DataProcessingPipe
    if use_gwaslab_rsids_convention:
        out_pipe = RenameColPipe(old_name="rsid", new_name="rsID")
    else:
        out_pipe = IdentityPipe()
    join_with_rsid_task = JoinDataFramesTask.create_from_result_df(
        asset_id=base_name + "_assign_rsids_via_dbsnp150",
        result_df_task=harmonize_task,
        reference_df_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_DIRECT_DOWNLOAD,
        left_on=["CHR", "POS", "EA", "NEA"],
        right_on=["int_chrom", "POS", "ALT", "REF"],
        out_format=ParquetOutFormat(),
        how="inner",
        df_1_pipe=CompositePipe(
            [
                CastPipe(target_column="EA", type=narwhals.dtypes.String(), new_col_name="EA"),
                CastPipe(target_column="NEA", type=narwhals.dtypes.String(), new_col_name="NEA"),
            ]
        ),
        backend="ibis",
        out_pipe=out_pipe,
    )
    return RSIDAssignmentTaskGroup(
        pre_harmonization_table_task=pre_harmonization_table_task,
        harmonize_task=harmonize_task,
        join_task=join_with_rsid_task,
    )
```

- [ ] **Step 2: Rewrite the Liu IBD harmonized asset**

Replace the whole content of liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized.py with:

```python
"""
Genome-reference harmonization of the Liu et al. IBD sumstats (lifted over to build 37),
against the UCSC hg19 genome and the 1000 Genomes EUR panel.

- Orient every variant so that NEA is the plus-strand reference allele.
- Resolve palindromic SNVs and ambiguous indels by panel frequency when the table's
  orientation is not fully reference-consistent, and drop them when the evidence is not
  decisive.
"""

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats import (
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)

LIU_ET_AL_2023_IBD_EUR_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
        asset_id="liu_et_al_2023_ibd_eur_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

LIU_ET_AL_2023_IBD_EUR_HARMONIZE = GenomeReferenceHarmonizationTask.create(
    asset_id="liu_et_al_2023_ibd_eur_genome_reference_harmonized",
    sumstats_task=LIU_ET_AL_2023_IBD_EUR_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
```

- [ ] **Step 3: Delete the dump asset and update its two consumers**

Run: `git rm mecfs_bio/assets/gwas/inflammatory_bowel_disease/liu_et_al_2023/processed_gwas_data/liu_et_al_2023_eur_37_harmonized_dump_to_parquet.py`

In both liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar.py and
liu_et_al_2023_eur_37_harmonized_assign_rsid_via_snp150_annovar_with_dups.py:

Replace

```python
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_37_harmonized_dump_to_parquet import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,
)
```

with

```python
from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats_harmonized import (
    LIU_ET_AL_2023_IBD_EUR_HARMONIZE,
)
```

and replace `result_df_task=LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET,` with
`result_df_task=LIU_ET_AL_2023_IBD_EUR_HARMONIZE,`.

Then run `grep -rn "LIU_ET_AL_2023_IBD_EUR_HARMONIZE_PARQUET\|harmonized_dump_to_parquet" mecfs_bio test_mecfs_bio docs`.

Expected: no Python matches. Documentation mentions, if any, get updated to the
new asset name.

- [ ] **Step 4: Rewrite the system test**

Replace test_mecfs_bio/system/test_harmonize_drop_ambiguous.py with:

```python
"""
System test of genome-reference harmonization on the first rows of DecodeME (lifted over
to build 37), using the real hg19 FASTA and 1000 Genomes EUR panel assets.
"""

import tempfile
from pathlib import Path

import polars as pl

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.filtered_snps_gwas_1 import (
    DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.head_pipe import HeadPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

_decode_me_first_rows = PipeDataFrameTask.create(
    source_task=DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
    asset_id="testing_decode_me_slice_task",
    out_format=ParquetOutFormat(),
    pipes=[HeadPipe(num_rows=10_000)],
)

_decode_me_first_rows_liftover_to_37 = GWASLabCreateSumstatsTask(
    df_source_task=_decode_me_first_rows,
    target_asset_id=AssetId("testing_first_rows_decode_me_gwas_1_sumstats_liftover_to_37"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
)

_pre_harmonization_table = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=_decode_me_first_rows_liftover_to_37,
    asset_id="testing_first_rows_pre_harmonization_table",
    sub_dir="processed",
)

_harmonized_task = GenomeReferenceHarmonizationTask.create(
    asset_id="testing_first_rows_genome_reference_harmonized",
    sumstats_task=_pre_harmonization_table,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)


def test_genome_reference_harmonization_orients_every_kept_variant_to_the_reference():
    with tempfile.TemporaryDirectory() as tempdirname:
        tempdir = Path(tempdirname)
        asset_root = tempdir / "asset_store"
        asset_root.mkdir(parents=True, exist_ok=True)
        runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=tempdir / "info_store.yaml",
            asset_root=asset_root,
        )
        assets = runner.run([_harmonized_task, UCSC_HG19_INDEXED_FASTA])
        harmonized = assets[_harmonized_task.asset_id]
        fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
        assert isinstance(harmonized, FileAsset)
        assert isinstance(fasta_asset, DirectoryAsset)
        table = pl.read_parquet(harmonized.path)
        assert table.height > 0
        fasta = IndexedFasta.open(fasta_asset.path)
        for chrom, group in table.group_by(GWASLAB_CHROM_COL):
            matches = reference_matches(
                fasta,
                chrom=int(chrom[0]),
                positions=group[GWASLAB_POS_COL].to_numpy(),
                alleles=group[GWASLAB_NON_EFFECT_ALLELE_COL],
            )
            assert matches.all()
```

- [ ] **Step 5: Enable disk cleanup in the system test workflow**

The run now holds the uncompressed hg19 FASTA (3.2 GB) plus the temporary panel
VCF and its sites TSV. In .github/workflows/system_test_harmonize.yml, change
`free-disk-space: "false"` to `free-disk-space: "true"`.

- [ ] **Step 6: Run green**

Run:

```bash
pixi r invoke green > /tmp/claude-green.log 2>&1; echo "EXIT=$?"
grep -aE "passed|failed|error|All checks passed" /tmp/claude-green.log | tail
```

Expected: EXIT=0 and no test failures. Import-time construction of every asset
module confirms the rewired graph.

- [ ] **Step 7: Run the system test locally**

Run: `pixi r python -m pytest test_mecfs_bio/system/test_harmonize_drop_ambiguous.py -q`

Expected: 1 passed. It builds reference assets in a temporary store, which
takes a while.

- [ ] **Step 8: Commit**

```bash
git add -A mecfs_bio/asset_generator/annovar_37_basic_rsid_assignment.py mecfs_bio/assets/gwas/inflammatory_bowel_disease test_mecfs_bio/system/test_harmonize_drop_ambiguous.py .github/workflows/system_test_harmonize.yml
git commit -m "Switch rsID-assignment chains to genome-reference harmonization"
```

### Task 12: Remove dead gwaslab harmonization code, survey trust (V5), final check

**Files:**
- Modify: mecfs_bio/build_system/task/gwaslab/gwaslab_create_sumstats_task.py
- Modify: mecfs_bio/asset_generator/ukbb_ppp_single_gwas_generator.py and mecfs_bio/assets/reference_data/ukbb_ppp_sumstats/rabgap1l/processed/ukbb_rabgap1l_sumstats_37.py (drop harmonize_options=None)
- Create: experiments/claude/genome_reference_harmonization/survey_trust.py
- Modify: experiments/claude/genome_reference_harmonization/validation_report.md

**Interfaces:**
- Consumes: RSIDAssignmentTaskGroup.pre_harmonization_table_task (Task 11), count_trust_evidence_genome_wide, decide_trust.

- [ ] **Step 1: Delete the dead code**

In gwaslab_create_sumstats_task.py delete:
- class HarmonizationOptions;
- class GWASLabVCFRef;
- function _do_harmonization;
- function _prune_unused_allele_categories;
- the harmonize_options field of GWASLabCreateSumstatsTask and of GwasLabTransformSpec;
- harmonize_options=self.harmonize_options in GWASLabCreateSumstatsTask.execute;
- the `if spec.harmonize_options is not None:` block in transform_gwaslab_sumstats;
- imports left unused by these deletions (keep gwaslab_download_ref_if_missing in gwaslab_util.py; the region-plot tasks use it).

Remove the `harmonize_options=None,` arguments in ukbb_ppp_single_gwas_generator.py
and ukbb_rabgap1l_sumstats_37.py.

Run: `grep -rn "HarmonizationOptions\|GWASLabVCFRef\b\|harmonize_options\|_do_harmonization\|_prune_unused_allele_categories" mecfs_bio test_mecfs_bio docs`

Expected: no matches. Documentation mentions get updated to describe
genome-reference harmonization.

- [ ] **Step 2: Write the trust survey script**

```python
"""
V5: trust decision of every rsID-assignment chain's pre-harmonization table.

Imports every asset module, collects RSIDAssignmentTaskGroup instances, and for each
group whose upstream gwaslab Sumstats pickle is already materialized, builds the
pre-harmonization table and reports trust counts. Groups whose inputs are not built are
listed as skipped, so nothing heavy is rebuilt.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.survey_trust \
    2>&1 | tee experiments/claude/genome_reference_harmonization/survey_trust.log
"""

import importlib
import pkgutil

import mecfs_bio.assets
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.annovar_37_basic_rsid_assignment import RSIDAssignmentTaskGroup
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import IndexedFasta
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    chromosomes_to_harmonize,
    count_trust_evidence_genome_wide,
    scan_sumstats_as_polars,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.trust import decide_trust
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

OPTIONS = GenomeReferenceHarmonizationOptions()


def _groups() -> list[RSIDAssignmentTaskGroup]:
    found: dict[str, RSIDAssignmentTaskGroup] = {}
    for module_info in pkgutil.walk_packages(mecfs_bio.assets.__path__, "mecfs_bio.assets."):
        module = importlib.import_module(module_info.name)
        for value in vars(module).values():
            if isinstance(value, RSIDAssignmentTaskGroup):
                found[value.harmonize_task.asset_id] = value
    return list(found.values())


def main() -> None:
    fasta_asset = DEFAULT_RUNNER.run([UCSC_HG19_INDEXED_FASTA])[UCSC_HG19_INDEXED_FASTA.asset_id]
    assert isinstance(fasta_asset, DirectoryAsset)
    fasta = IndexedFasta.open(fasta_asset.path)
    for group in _groups():
        table_task = group.pre_harmonization_table_task
        assert isinstance(table_task, GwasLabSumstatsToTableTask)
        pickle_path = DEFAULT_RUNNER.meta_to_path(table_task.source_sumstats_task.meta)
        if not pickle_path.exists():
            print(f"SKIPPED (input not materialized): {group.harmonize_task.asset_id}")
            continue
        table_asset = DEFAULT_RUNNER.run([table_task])[table_task.asset_id]
        sumstats = scan_sumstats_as_polars(table_asset, table_task.meta, IdentityPipe())
        chromosomes = chromosomes_to_harmonize(sumstats, fasta, OPTIONS)
        counts = count_trust_evidence_genome_wide(sumstats, chromosomes, fasta, OPTIONS)
        print(
            f"{group.harmonize_task.asset_id}: trusted={decide_trust(counts, OPTIONS)} {counts}"
        )


if __name__ == "__main__":
    main()
```

- [ ] **Step 3: Run the survey and record it**

Run: `pixi r python -m experiments.claude.genome_reference_harmonization.survey_trust 2>&1 | tee experiments/claude/genome_reference_harmonization/survey_trust.log`

Expected: one line per group, either a trust line or SKIPPED. Copy the table
into validation_report.md.

- [ ] **Step 4: Final green**

Run:

```bash
pixi r invoke green > /tmp/claude-green.log 2>&1; echo "EXIT=$?"
grep -aE "passed|failed|error|All checks passed" /tmp/claude-green.log | tail
```

Expected: EXIT=0, no failures.

- [ ] **Step 5: Commit**

```bash
git add -A mecfs_bio test_mecfs_bio experiments/claude/genome_reference_harmonization
git commit -m "Remove gwaslab harmonization code and survey trust decisions (V5)"
```

- [ ] **Step 6: Prepare the PR description (only when the user asks for a PR)**

Include:
- the goal;
- the clarifications to the spec (Global Constraints);
- the V1-V5 summary from validation_report.md;
- the analyses whose inputs change: every chain built by annovar_37_basic_rsid_assignment, and Liu IBD. Indel sets and palindrome handling shift, and asset ids change, so downstream assets rebuild;
- the out-of-scope follow-ups from the spec.

---

## Self-Review Notes

**Spec coverage**

| Spec requirement | Task |
|---|---|
| Pinned FASTA and panel assets, with DiscardDepsWrapper | Tasks 1-3 |
| Duplicate-key handling in the panel | Task 2 |
| Allele classes and reverse-strand rules | Task 4 |
| Palindromic MNP (trusted kept, untrusted dropped) | Task 4 |
| Trust (100%, minimum counts) | Tasks 4 and 6 |
| Palindrome rules | Task 5 |
| Stringent indel rules and symmetry | Task 6 |
| Flip registry (A1FREQ, bounds, unregistered columns, extra rules) | Task 4 |
| Uniqueness assertion | Task 4 |
| Streaming per chromosome; backend assertion | Task 4 |
| Excluded MT | Task 4 |
| V1/V2 parity | Task 8 |
| V3 tuning | Task 9 |
| V4 memory | Task 10 |
| V5 trust survey | Task 12 |
| Rollout (call sites, system test) | Task 11 |
| Rollout (dead code) | Task 12 |

Deviations are listed under Global Constraints, "Clarifications to the spec".

**Type consistency.** These names are used identically across tasks:
- decide_palindrome_strands returns pl.Series;
- decide_ambiguous_indels returns a two-column pl.DataFrame;
- ParquetPanelLoader(panel_path, chrom) is called with positions;
- resolve_chromosome_rows(sumstats, context, panel_path);
- RSIDAssignmentTaskGroup.pre_harmonization_table_task.
