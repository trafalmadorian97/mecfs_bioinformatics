"""
Pinned constants for genome-wide fine-mapping (GWFM) with the GCTB binary's SBayesRC
model (gctb --gwfm RC).
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


def gwfm_reference_prefix(version: str) -> str:
    """Return the bucket-relative S3 key prefix for a staged GWFM reference.
    The trailing slash lets a filename be appended directly .
    """
    return f"sbayesrc/reference/{version}/"


@frozen
class GCTBReferenceBundleFile:
    filename: str
    source_url: str
    size_bytes: int
    sha256: str | None  # None until first staging records it, then committed


# Inner names of the pinned reference files
# Contents of zips:
# ukbEUR_13M_FullLDM.zip -> ldm13M/ (snp.info, ldm.info, rsq0.5.pwld, block*.ldm.bin);
# annot_baseline2.2_13M.zip -> the single file annot_baseline2.2_13M.txt.
GWFM_LDM_ZIP_NAME: str = "ukbEUR_13M_FullLDM.zip"
GWFM_ANNOT_ZIP_NAME: str = "annot_baseline2.2_13M.zip"
GWFM_LDM_DIR_NAME: str = "ldm13M"
GWFM_ANNOT_FILE_NAME: str = "annot_baseline2.2_13M.txt"
GWFM_GENE_MAP_FILE_NAME: str = "gene_map_hg38_hg19.txt"
GWFM_PWLD_RELPATH: str = "ldm13M/rsq0.5.pwld"

GWFM_REFERENCE_BUNDLE: tuple[GCTBReferenceBundleFile, ...] = (
    GCTBReferenceBundleFile(
        GWFM_LDM_ZIP_NAME,
        "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/ukbEUR_13M_FullLDM.zip",
        206566549726,
        None,
    ),
    GCTBReferenceBundleFile(
        "ref_b37_1588blocks.pos",
        "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/LD/Imputed13M/ref_b37_1588blocks.pos",
        40058,
        None,
    ),
    GCTBReferenceBundleFile(
        GWFM_ANNOT_ZIP_NAME,
        "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/GWFM/Annotation/annot_baseline2.2_13M.zip",
        557227563,
        None,
    ),
    GCTBReferenceBundleFile(
        GWFM_GENE_MAP_FILE_NAME,
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


# Remote container layout: the /work mount holds staged reference files under work/ref,
# the .ma sumstats at work/sumstats.ma, and all gctb outputs under work/out.
REMOTE_REF_DIR: str = "work/ref"
REMOTE_OUT_DIR: str = "work/out"
REMOTE_MA_PATH: str = "work/sumstats.ma"
MATCHED_LDM_OUT: str = "work/out/matched_ldm"
GWFM_OUT_PREFIX: str = "work/out/gwfm"

# Credible-set posterior-inclusion / posterior-enrichment probability thresholds
DEFAULT_PIP: float = 0.9
DEFAULT_PEP: float = 0.7

MARKER_VERSION_KEY: str = "version"
MARKER_PREFIX_KEY: str = "s3_prefix"
MARKER_FILES_KEY: str = "files"

# Column order for a GCTB/COJO .ma summary-statistics file (tab-separated, with
# header). SNP carries the variant id (rsID) and must match the variant ids used by
# the LD reference (Task 1's GWFM_REFERENCE_BUNDLE), since GCTB joins the .ma to the
# LD matrix by this column.
COJO_MA_COLUMNS: tuple[str, ...] = ("SNP", "A1", "A2", "freq", "b", "se", "p", "N")
