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
