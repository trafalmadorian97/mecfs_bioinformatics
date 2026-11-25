"""Column names of the docs-facing SUSIE-PolyFun explanation tables.

Kept free of heavy imports so the docs macros (docs_macros) can import these names
to describe the columns in figure captions.
"""

DISP_CHR = "chr"
DISP_POS = "pos"
DISP_EA = "ea"
DISP_NEA = "nea"
DISP_CS_PF = "cs_pf"
DISP_CS_U = "cs_u"
DISP_PIP_PF = "pip_pf"
DISP_PIP_U = "pip_u"
DISP_LIFT = "lift"
# Prefix on the detailed table's per-family contrast columns, e.g. annot_coding.
DISP_ANNOT_PREFIX = "annot_"

# Per-annotation context columns on the per-variant table: the annotation
# coefficient gamma_raw_c (see ridge_weights_task), and abar_c (the uniform-run
# PIP-weighted mean of the annotation over the locus variants).
DISP_GAMMA = "gamma"
DISP_ALPHA_BAR = "alpha_bar"

# Columns of the prior comparison table, which sets an external and an internal
# PolyFun prior's SUSIE runs side by side with the uniform-prior run.
CS_PF_EXT_COL = "cs_pf_ext"
CS_PF_INT_COL = "cs_pf_int"
CS_U_COL = "cs_u"
PIP_PF_EXT_COL = "pip_pf_ext"
PIP_PF_INT_COL = "pip_pf_int"
PIP_U_COL = "pip_u"
LIFT_EXT_COL = "lift_ext"
LIFT_INT_COL = "lift_int"


def secondary_pos_display_col(build_label: str) -> str:
    """Name of the display column holding a secondary-build position, e.g. pos_hg38."""
    return f"{DISP_POS}_{build_label}"
