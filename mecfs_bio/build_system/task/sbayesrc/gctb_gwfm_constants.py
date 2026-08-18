"""
Pinned constants for genome-wide fine-mapping (GWFM) with the GCTB binary's SBayesRC
model (gctb --gwfm RC).
"""

from attrs import frozen

# Upstream base version the fork descends from (informational; the image is now
# built from source, not from the pinned upstream binary).
GCTB_VERSION: str = "2.5.5"

# The gctb Docker image is built FROM SOURCE against our GCTB fork so we can
# patch upstream bugs (the make-ldm-eigen concurrency race, etc.); see
# plan_for_dummy_run_revised.txt. Bump GCTB_FORK_REF after each fork commit --
# GCTB_IMAGE_TAG derives from it, so a new binary always gets a new, unambiguous
# tag that `docker pull` fetches fresh.
GCTB_FORK_REPO: str = "https://github.com/trafalmadorian97/GCTB_traf_fork.git"
GCTB_FORK_REF: str = "86d3a03ac49dab5b18fc3fce42dd933e138b1db0"  # branch fix/gwfm-segfaults
GCTB_IMAGE_TAG: str = f"{GCTB_VERSION}-fork-{GCTB_FORK_REF[:12]}"


def gctb_image_build_args() -> list[str]:
    """
    Docker --build-arg pairs that build the gctb image from the pinned fork source.

    Shared by every docker-build call site (the publish task, the toy-reference
    builder, and the system tests) so they stay in lockstep on the fork/ref pin.
    """
    return [
        "--build-arg",
        f"GCTB_FORK_REPO={GCTB_FORK_REPO}",
        "--build-arg",
        f"GCTB_FORK_REF={GCTB_FORK_REF}",
    ]


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
# The 13M reference unzips to ~260 GiB; peak disk is that plus the make-ldm-eigen
# output (whose size is not published). Even with the zip deleted right after unzip
# (see _build_commands), 500 GiB was insufficient. 1000 GiB is generous headroom;
# EBS is a few cents for a several-hour run, so over-provisioning is cheap insurance.
DEFAULT_DISK_GB: int = 1000


# Remote container layout: the /work mount holds staged reference files under work/ref,
# the .ma sumstats at work/sumstats.ma, and the small gctb result files under work/out
# (the only directory retrieved back). The resized/eigen-decomposed LD matrices are a
# large intermediate (tens of GiB) consumed only by the gwfm step on the same instance,
# so they are written to work/matched_ldm -- a SIBLING of work/out, deliberately NOT
# under it -- to keep them off the output round-trip (instance -> scratch S3 -> local).
REMOTE_REF_DIR: str = "work/ref"
REMOTE_OUT_DIR: str = "work/out"
REMOTE_MA_PATH: str = "work/sumstats.ma"
MATCHED_LDM_OUT: str = "work/matched_ldm"
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
