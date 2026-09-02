"""Column names for the parsed Eagle genetic-map parquet.

The parsed map is keyed by (GWASLAB_CHROM_COL, GMAP_POS_COL) and carries the
per-position recombination rate (cM/Mb) and the cumulative genetic-map position
(cM). These names are shared by the task that builds the parquet and the
explainability plot's recombination track that reads it.
"""

GMAP_POS_COL = "POS"
GMAP_RATE_COL = "recomb_rate_cm_per_mb"
GMAP_CM_COL = "genetic_map_cm"
