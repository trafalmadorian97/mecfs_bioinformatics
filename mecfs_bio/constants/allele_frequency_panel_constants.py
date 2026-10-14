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

# gnomAD genetic-ancestry group codes, with gnomAD's own names. Source: the gnomAD FAQ
# "How are genetic ancestry group names abbreviated?",
# https://gnomad.broadinstitute.org/help/how-are-genetic-ancestry-group-names-abbreviated
# (text in https://github.com/broadinstitute/gnomad-browser/blob/main/browser/help/faq/technical-details/how-are-genetic-ancestry-group-names-abbreviated.md).
#   afr: African/African American
#   ami: Amish
#   amr: Admixed American
#   asj: Ashkenazi Jewish
#   eas: East Asian
#   fin: European (Finnish)
#   mid: Middle Eastern
#   nfe: European (non-Finnish)
#   nfe_est: Estonian (nfe subgroup, v2)
#   nfe_nwe: North-western European (nfe subgroup, v2), e.g. UK and Nordic cohorts
#   nfe_onf: Other non-Finnish European (nfe subgroup, v2)
#   nfe_seu: Southern European (nfe subgroup, v2)
#   oth: Other (v2)
#   remaining: Remaining individuals (v4; the successor of oth)
#   sas: South Asian
# Groups present in gnomAD v2.1.1 genomes (the hg19 panel).
GnomadV2Group = Literal[
    "afr",
    "amr",
    "asj",
    "eas",
    "fin",
    "nfe",
    "nfe_est",
    "nfe_nwe",
    "nfe_onf",
    "nfe_seu",
    "oth",
]
# Groups present in gnomAD v4.1 genomes (the hg38 panel).
GnomadV4Group = Literal[
    "afr",
    "ami",
    "amr",
    "asj",
    "eas",
    "fin",
    "mid",
    "nfe",
    "remaining",
    "sas",
]
GnomadGroup = GnomadV2Group | GnomadV4Group
PanUkbbGroup = Literal["ukb_afr", "ukb_amr", "ukb_csa", "ukb_eas", "ukb_eur", "ukb_mid"]
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
