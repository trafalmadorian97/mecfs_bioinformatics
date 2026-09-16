# Delta: suspicious-indel trust check

Amends 2026-09-14-genome-reference-harmonization-{design,plan}.md. Folded into
both on 2026-09-16 (design: trust-decision gate, no-EAF rule, monomorphic-record
rule, revised V3, decisions 6-9; plan: options fields, Task 6 step 4b, revised
Task 9). This document is retained as the detailed rationale and full code.
Motivation and evidence: v3_trust_investigation_summary.md.

## Why

Passing the 100% SNV+indel consistency check does not guarantee that a table's
ambiguous indels (class BOTH: both orientations match the FASTA) follow the
reference convention. The NEA==REF test is blind to orientation there. Measured
share of ambiguous indels whose panel frequency decisively contradicts the source
orientation (distance 0.02, margin 0.3):

| Dataset | rate | want |
|---|---|---|
| DecodeME build 38 (raw, imputed to a GRCh38 WGS panel) | 0.0014% | trusted |
| MVP myocardial infarction (GWAS Catalog harmonised) | 0.026% | untrusted |
| Bellenguez Alzheimer's (harmonised) | 2.2% | untrusted |
| Kerrebijn fibromyalgia (harmonised) | 8.3% | untrusted |

A ~20x gap separates DecodeME from the harmonised files. So trust should also
require the proportion of "suspicious" ambiguous indels to be low.

## Definition

A **suspicious** ambiguous indel is one that `decide_ambiguous_indels` resolves
to ACTION_SWAP: its EAF matches the complement of the opposite-orientation panel
AF decisively (the flip reading fits within indel_max_af_distance and beats the
keep reading by indel_min_af_margin). This is exactly the user's "EAF matches the
complement of the opposite orientation's alt frequency", reusing existing code.

A **checkable** ambiguous indel is one with EAF present and at least one panel
record (keep or flip) -- i.e. one `decide_ambiguous_indels` does not drop for
NO_EAF or NOT_IN_PANEL. suspicious is a subset of checkable.

**suspicious fraction** = suspicious / checkable.

The suspicious-fraction test is a positive misorientation detector over tables
that carry EAF: it can only fail trust, and a too-sparse checkable set does not
fail it. Separately, a table with no EAF column cannot assess ambiguous-indel
orientation at all, so it is refused trust outright (decision, 2026-09-16); its
ambiguous indels are then dropped on the untrusted path (NO_EAF).

## Options (options.py)

Add two fields:

```python
    max_suspicious_indel_fraction: float = 1e-4
    min_checkable_ambiguous_indels: int = 100
```

Docstring lines:

    max_suspicious_indel_fraction: a table is trusted only if at most this
        fraction of its checkable ambiguous indels (both orientations on the
        FASTA, with EAF and a panel record) are suspicious -- their EAF matches
        the complement of the opposite orientation's panel frequency, evidence
        the source is not perfectly reference-aligned. The suspicion test uses
        the same indel_max_af_distance and indel_min_af_margin.
    min_checkable_ambiguous_indels: the suspicious-fraction test is applied only
        when at least this many checkable ambiguous indels exist; below it the
        signal is too sparse to judge and does not affect trust.

post_init asserts:

```python
        assert 0 < self.max_suspicious_indel_fraction < 1
        assert self.min_checkable_ambiguous_indels >= 1
```

Default value note: 1e-4 sits in the gap (DecodeME ~1.4e-5, next dataset
~2.6e-4). Confirmed/adjusted by the revised V3.

## trust.py

Add a counts container and a counter, and extend decide_trust.

```python
from mecfs_bio.build_system.task.genome_reference_harmonization.allele_classes import (
    ALLELE_CLASS_COL,
    CLASS_INDEL_BOTH,
    ...
)
from mecfs_bio.build_system.task.genome_reference_harmonization.ambiguous_indels import (
    INDEL_ACTION_COL,
    INDEL_DROP_REASON_COL,
    decide_ambiguous_indels,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.outcomes import (
    ACTION_SWAP,
    DROP_AMBIGUOUS_INDEL_AF_INDECISIVE,
    DROP_AMBIGUOUS_INDEL_AF_MISMATCH,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

_SUSPICIOUS = "suspicious"
_CHECKABLE = "checkable"


@frozen
class SuspiciousIndelCounts:
    suspicious: int  # ambiguous indels whose panel frequency contradicts the source orientation
    checkable: int  # ambiguous indels with EAF and at least one panel record

    def __add__(self, other: "SuspiciousIndelCounts") -> "SuspiciousIndelCounts":
        return SuspiciousIndelCounts(
            suspicious=self.suspicious + other.suspicious,
            checkable=self.checkable + other.checkable,
        )

    @classmethod
    def zero(cls) -> "SuspiciousIndelCounts":
        return cls(suspicious=0, checkable=0)

    @property
    def fraction(self) -> float:
        return self.suspicious / self.checkable if self.checkable else 0.0


@frozen
class TrustEvidence:
    counts: TrustCounts
    suspicious: SuspiciousIndelCounts
    eaf_present: bool


def count_suspicious_indels(
    classified: pl.DataFrame,
    panel: pl.DataFrame,
    options: GenomeReferenceHarmonizationOptions,
) -> SuspiciousIndelCounts:
    """Suspicious and checkable ambiguous indels in a classify_alleles frame that carries EAF."""
    ambiguous = classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH)
    if ambiguous.height == 0 or GWASLAB_EFFECT_ALLELE_FREQ_COL not in ambiguous.columns:
        return SuspiciousIndelCounts.zero()
    rows = ambiguous.select(
        GWASLAB_POS_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_NON_EFFECT_ALLELE_COL,
        pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).cast(pl.Float64),
    )
    decisions = decide_ambiguous_indels(rows, panel, options)
    reason = pl.col(INDEL_DROP_REASON_COL)
    checkable = reason.is_null() | reason.is_in(
        [DROP_AMBIGUOUS_INDEL_AF_MISMATCH, DROP_AMBIGUOUS_INDEL_AF_INDECISIVE]
    )
    summary = decisions.select(
        (pl.col(INDEL_ACTION_COL) == ACTION_SWAP).sum().alias(_SUSPICIOUS),
        checkable.sum().alias(_CHECKABLE),
    ).row(0, named=True)
    return SuspiciousIndelCounts(
        suspicious=int(summary[_SUSPICIOUS]), checkable=int(summary[_CHECKABLE])
    )


def decide_trust(
    evidence: TrustEvidence, options: GenomeReferenceHarmonizationOptions
) -> bool:
    counts = evidence.counts
    consistent = (
        counts.inconsistent_snvs == 0
        and counts.inconsistent_indels == 0
        and counts.consistent_snvs >= options.min_checkable_snvs
        and counts.consistent_indels >= options.min_checkable_indels
    )
    if not consistent:
        return False
    if not evidence.eaf_present:
        # No EAF column: ambiguous-indel orientation cannot be assessed, so refuse
        # trust (decision, 2026-09-16). Untrusted resolution then drops them (NO_EAF).
        return False
    if evidence.suspicious.checkable < options.min_checkable_ambiguous_indels:
        return True
    return evidence.suspicious.fraction <= options.max_suspicious_indel_fraction
```

decide_trust now takes a TrustEvidence (breaking change to its two call sites:
the task module and survey_trust.py). count_suspicious_indels is unchanged: when
EAF is absent it returns zero(), and the EAF-absent refusal is driven by
evidence.eaf_present, not by the counts.

## Task module (genome_reference_harmonization_task.py)

Pass 1 now also reads EAF and consults the panel at BOTH-indel positions.

```python
def _ambiguous_positions(classified: pl.DataFrame) -> pl.Series:
    return classified.filter(pl.col(ALLELE_CLASS_COL) == CLASS_INDEL_BOTH)[GWASLAB_POS_COL]


def count_trust_evidence_genome_wide(
    sumstats: pl.LazyFrame,
    chromosomes: Sequence[int],
    fasta: IndexedFasta,
    panel_path: Path,
    options: GenomeReferenceHarmonizationOptions,
) -> TrustEvidence:
    names = sumstats.collect_schema().names()
    eaf_present = GWASLAB_EFFECT_ALLELE_FREQ_COL in names
    columns = [c for c in [*_KEY_COLUMNS, GWASLAB_EFFECT_ALLELE_FREQ_COL] if c in names]
    counts = TrustCounts.zero()
    suspicious = SuspiciousIndelCounts.zero()
    for chrom in chromosomes:
        frame = (
            sumstats.filter(pl.col(GWASLAB_CHROM_COL) == chrom)
            .select(columns)
            .collect(engine="streaming")
        )
        classified = classify_alleles(
            prepare_alleles(frame).filter(valid_alleles_expr()),
            fasta=fasta,
            chrom=chrom,
            max_gather_bytes=options.max_gather_bytes,
        )
        counts = counts + count_trust_evidence(classified)
        if eaf_present:
            panel = ParquetPanelLoader(panel_path=panel_path, chrom=chrom)(
                _ambiguous_positions(classified)
            )
            suspicious = suspicious + count_suspicious_indels(classified, panel, options)
    return TrustEvidence(counts=counts, suspicious=suspicious, eaf_present=eaf_present)
```

execute():

```python
        evidence = count_trust_evidence_genome_wide(
            sumstats, chromosomes, panel_asset.path, self.options
        )
        trusted = decide_trust(evidence, self.options)
        logger.info(
            "genome-reference harmonization trust decision",
            trusted=trusted,
            eaf_present=evidence.eaf_present,
            counts=attrs.asdict(evidence.counts),
            suspicious=attrs.asdict(evidence.suspicious),
        )
```

Note the argument order change (panel_asset.path added before options).
ParquetPanelLoader must be defined above count_trust_evidence_genome_wide, or the
positions helper hoisted; a small reorder.

## Tests (test_genome_reference_harmonization_task.py)

TEST_OPTIONS in the fixtures gains min_checkable_ambiguous_indels=1 so a single
row can drive the fraction test.

```python
def test_suspicious_ambiguous_indels_make_a_table_untrusted(tmp_path: Path) -> None:
    # pos 5 is T at the head of a G run, so T/TG is class BOTH. The panel has only the
    # opposite reading (REF=T, ALT=TG) at AF 0.9, and EAF 0.1 = 1 - 0.9, so the row is
    # suspicious. With min_checkable_ambiguous_indels=1 that makes the table untrusted,
    # and the ambiguous indel is swapped rather than kept in source orientation.
    ambiguous = Variant(pos=5, ea="T", nea="TG", eaf=0.1)
    panel = [PanelRecord(pos=5, ref="T", alt="TG", af=0.9)]
    result = run_harmonization(
        tmp_path / "run",
        sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, ambiguous]),
        panel=panel,
    )
    swapped = row_at(result, 5)
    assert (swapped[EA], swapped[NEA]) == ("TG", "T")
```

```python
def test_table_without_eaf_is_untrusted(tmp_path: Path) -> None:
    # With no EAF column, ambiguous-indel orientation cannot be assessed, so the table
    # is refused trust and its ambiguous indel is dropped (NO_EAF) rather than kept.
    ambiguous = Variant(pos=5, ea="T", nea="TG")
    frame = sumstats_frame([CONSISTENT_SNV, CONSISTENT_INDEL, ambiguous]).drop(
        GWASLAB_EFFECT_ALLELE_FREQ_COL
    )
    result = run_harmonization(tmp_path / "run", frame)
    assert 5 not in positions(result)
    assert set(positions(result)) >= {1, 21}
```

The existing test_trusted_table_keeps_ambiguous_indel_in_source_orientation stays
green: with no panel the sole ambiguous indel is NOT_IN_PANEL, so checkable = 0,
the sparse branch applies, and the table stays trusted and keeps (T, TG).

## Ambiguous-indel rule change (Task 6, decision 2026-09-16)

Two decisions on decide_ambiguous_indels, which both the untrusted resolution and
the suspicion counter use:

- **Drop monomorphic panel records.** Filter the panel to strictly polymorphic
  records at the top of decide_ambiguous_indels, so a record that predicts EAF 0
  or 1 is treated as absent:

  ```python
      panel = panel.filter((pl.col(PANEL_AF_COL) > 0) & (pl.col(PANEL_AF_COL) < 1))
  ```

  On build-38 DecodeME this drops the AF-0 group of contradicting swaps (11 -> 5).
  It changes suspicious counts, so the max_suspicious_indel_fraction calibration
  in V3 uses the filtered panel.

- **Do not require both orientations.** Keep the current single-record logic: a
  reading is chosen when one record fits (and, if both exist, beats the other by
  the margin). The require-both variant was rejected because it does not improve
  safety on mis-oriented tables (Kerrebijn still ~100k wrong) and only lowers
  yield.

Note: some existing Task-6 ambiguous-indel tests use AF-0 records to represent a
misfitting record; with the filter those become absent (NOT_IN_PANEL) rather than
misfit (AF_MISMATCH). Revisit those fixtures when folding this in. The palindrome
rule is out of scope and keeps its current panel handling.

## Ordering in the plan

count_suspicious_indels depends on decide_ambiguous_indels, which Task 6 creates.
So implement this delta as a new step at the end of Task 6 (after ambiguous
indels exist), not in Task 4. Task 4 ships decide_trust(counts, options) as
today; Task 6 extends it to decide_trust(counts, suspicious, options), adds the
two options, count_suspicious_indels/SuspiciousIndelCounts/TrustEvidence, rewires
count_trust_evidence_genome_wide and execute, and adds the test above.

## Revised V3 (Task 9)

The original V3 goal (a zero-wrong grid cell) is unreachable: no
(distance, margin) gives zero wrong on build-38 DecodeME, because a trusted table
still carries ~11 ambiguous indels that differ only in panel representation
(UKB WGS vs 1000 Genomes) in repeats. Revise V3 to:

1. Record that finding and the cross-dataset investigation
   (v3_trust_investigation_summary.md, v3_other_truth_sets.py).
2. Set the untrusted-path resolution defaults (indel_max_af_distance,
   indel_min_af_margin) to the most conservative sensible cell; recommend
   0.02 and 0.3 (fewest wrong swaps in the grid). Untrusted tables are exactly
   those whose orientation is unreliable, so aggressive dropping is appropriate.
3. Calibrate max_suspicious_indel_fraction so DecodeME build 38 is trusted and
   the GWAS Catalog harmonised datasets (MVP, Bellenguez, Kerrebijn) are
   untrusted; record each dataset's suspicious fraction and confirm 1e-4.

## Design-spec text edits

- Trust decision (pass 1): add a third gate after checkable SNVs and indels:
  "Suspicious ambiguous indels: a BOTH indel is suspicious when its EAF matches
  the complement of the opposite orientation's panel AF decisively (ACTION_SWAP).
  trusted additionally requires the suspicious fraction over checkable ambiguous
  indels to be <= max_suspicious_indel_fraction, tested only when at least
  min_checkable_ambiguous_indels are checkable." Note that the check is a
  positive detector: too-sparse or EAF-absent sets never fail trust.
- Pass 1 streaming section: pass 1 now collects CHR, POS, EA, NEA and EAF, and
  reads the panel at BOTH-indel positions to run the suspicion test; peak memory
  is still one chromosome. Pass 2 still reads the panel only when untrusted.
- Background: add a one-line pointer to the investigation and that DecodeME's
  build-38 alleles come from a GRCh38 UKB WGS imputation panel (supplement).
- Replace the V3 description with the revised V3 above.

## Open decisions

1. No-EAF tables (e.g. PGC schizophrenia): RESOLVED (2026-09-16) -- refuse trust
   when the EAF column is absent, via evidence.eaf_present in decide_trust. Such a
   table takes the untrusted path, where its ambiguous indels are dropped (NO_EAF).
2. Default max_suspicious_indel_fraction = 1e-4 and min_checkable_ambiguous_indels
   = 100 are provisional; V3 confirms.
3. Ambiguous-indel rule: RESOLVED (2026-09-16) -- drop 0/1-AF panel records, and
   do not require both orientations. Folded in above (Ambiguous-indel rule change).
