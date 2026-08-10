"""
Pinned constants for genome-wide fine-mapping (GWFM) with the GCTB binary's SBayesRC
model (gctb --gwfm RC).

Values here were recorded during reconnaissance of the gctbhub site
(https://gctbhub.cloud.edu.au/software/gctb/) on 2026-08-09; see
experiments/claude/design_specs/2026-08-08-gwfm-recon-findings.md for the sourced facts
and reasoning behind each value, including why the GWFM_REFERENCE_BUNDLE excludes the
precomputed eigen/blk1588_eigen.tar.gz download.

sha256 fields on GWFM_REFERENCE_BUNDLE entries are None until the first real staging run
(a later task) records and commits the true checksums; only filename, source_url, and
size_bytes (from a HEAD request against gctbhub) are pinned now.
"""

from attrs import frozen

GCTB_VERSION: str = "2.5.5"
GCTB_BINARY_URL: str = (
    "https://gctbhub.cloud.edu.au/software/gctb/download/gctb_2.5.5_Linux.zip"
)
GCTB_BINARY_SHA256: str = (
    "ccc2752e1bc0d4a210bd9592d07c44468891677f8b2d5fc0a37aa2119c24c61f"
)

GWFM_REFERENCE_VERSION: str = "Imputed13M/v1"


@frozen
class ReferenceBundleFile:
    filename: str
    source_url: str
    size_bytes: int
    sha256: str | None  # None until first staging records it, then committed


GWFM_REFERENCE_BUNDLE: tuple[ReferenceBundleFile, ...] = (
    ReferenceBundleFile(
        "ukbEUR_13M_FullLDM.zip",
        "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/ukbEUR_13M_FullLDM.zip",
        206566549726,
        None,
    ),
    ReferenceBundleFile(
        "ref_b37_1588blocks.pos",
        "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/ref_b37_1588blocks.pos",
        40058,
        None,
    ),
    ReferenceBundleFile(
        "annot_baseline2.2_13M.zip",
        "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/Annotation/annot_baseline2.2_13M.zip",
        557227563,
        None,
    ),
    ReferenceBundleFile(
        "gene_map_hg38_hg19.txt",
        "https://gctbhub.cloud.edu.au/software/gctb/download/gene_map_hg38_hg19.txt",
        5155268,
        None,
    ),
)

# Fields filled at run time: {ldm} {ma} {annot} {gene_map} {out} {threads} {pip} {pep} {pwld}
GCTB_MAKE_LDM_EIGEN_TEMPLATE: str = "gctb --ldm {ldm} --gwas-summary {ma} --make-ldm-eigen --thread {threads} --out {out}"
GCTB_GWFM_TEMPLATE: str = (
    "gctb --gwfm RC --ldm-eigen {ldm} --gwas-summary {ma} --annot {annot} "
    "--gene-map {gene_map} --thread {threads} --out {out}"
)
GCTB_CS_TEMPLATE: str = (
    "gctb --cs --pwld-file {pwld} --pip {pip} --pep {pep} --gene-map {gene_map} "
    "--mcmc-samples {out} --out {out}"
)

DEFAULT_MEMORY_GB: int = 192
DEFAULT_VCPUS: int = 24
DEFAULT_DISK_GB: int = 500

MARKER_VERSION_KEY: str = "version"
MARKER_PREFIX_KEY: str = "s3_prefix"
MARKER_FILES_KEY: str = "files"

# Column order for a GCTB/COJO .ma summary-statistics file (tab-separated, with
# header). SNP carries the variant id (rsID) and must match the variant ids used by
# the LD reference (Task 1's GWFM_REFERENCE_BUNDLE), since GCTB joins the .ma to the
# LD matrix by this column.
COJO_MA_COLUMNS: tuple[str, ...] = ("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
