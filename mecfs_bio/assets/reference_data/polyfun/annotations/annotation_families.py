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
