"""Column names of the annotation-weights table (one row per baseline-LF annotation).

Written by RidgeAnnotationWeightsTask and L2SldscSnpvarTask. gamma_raw is the
annotation coefficient on the raw annotation scale; gamma_standardized is the same
coefficient per standard deviation of the annotation.
"""

ANNOTATION_COL = "annotation"
GAMMA_RAW_COL = "gamma_raw"
GAMMA_STANDARDIZED_COL = "gamma_standardized"
FAMILY_COL = "family"
