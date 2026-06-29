"""
Write GWAS summary statistics in COJO format, the input format consumed by
SBayesRC and polypwas.

The COJO format is a TAB delimited text table with the columns

    SNP A1 A2 freq b se p N

where SNP is the variant id used to match the LD reference, A1 is the effect
allele, A2 is the other allele, freq is the frequency of A1, b and se are the
marginal effect and its standard error on the original (unstandardized) scale,
p is the p value, and N is the per-SNP sample size.

This single writer feeds the SBayesRC ma_file as well as the polypwas --pqtl and
--gwas inputs, so all three stay in sync with the producer's column conventions.
"""

from pathlib import Path

import narwhals as nw

from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    PhenotypeInfo,
    QuantPhenotype,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_P_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL,
)

COJO_SNP_COL = "SNP"
COJO_A1_COL = "A1"
COJO_A2_COL = "A2"
COJO_FREQ_COL = "freq"
COJO_B_COL = "b"
COJO_SE_COL = "se"
COJO_P_COL = "p"
COJO_N_COL = "N"

COJO_COLUMN_ORDER = [
    COJO_SNP_COL,
    COJO_A1_COL,
    COJO_A2_COL,
    COJO_FREQ_COL,
    COJO_B_COL,
    COJO_SE_COL,
    COJO_P_COL,
    COJO_N_COL,
]


def _with_sample_size_column(
    frame: nw.LazyFrame, phenotype: PhenotypeInfo
) -> nw.LazyFrame:
    """
    Add the COJO N column from the phenotype info.

    For a quantitative phenotype with no recorded total sample size the incoming
    frame must already carry an N column (gwaslab names it N as well).  For a
    binary phenotype we use the effective sample size, consistent with how the
    rest of the codebase prepares per-SNP sample sizes.
    """
    if isinstance(phenotype, QuantPhenotype):
        if phenotype.total_sample_size is not None:
            return frame.with_columns(
                nw.lit(phenotype.total_sample_size).alias(COJO_N_COL)
            )
        # Fail fast (shift-left) rather than silently produce a frame missing N.
        columns = frame.collect_schema().names()
        assert COJO_N_COL in columns, (
            f"QuantPhenotype has no total_sample_size and the input frame has no "
            f"'{COJO_N_COL}' column (columns: {columns}); cannot determine per-SNP N."
        )
        return frame
    if isinstance(phenotype, BinaryPhenotypeSampleInfo):
        return frame.with_columns(
            nw.lit(phenotype.effective_sample_size).alias(COJO_N_COL)
        )
    raise ValueError(f"Unknown phenotype info: {phenotype}")


def write_cojo_ma(
    frame: nw.LazyFrame, phenotype: PhenotypeInfo, out_path: Path
) -> Path:
    """
    Write a gwaslab-format lazy frame to a COJO format .ma file at out_path.

    The frame must contain the gwaslab columns rsID, EA, NEA, EAF, BETA, SE and P.
    The N column is filled from the phenotype info (see _with_sample_size_column).
    """
    renamed = frame.with_columns(
        nw.col(GWASLAB_RSID_COL).alias(COJO_SNP_COL),
        nw.col(GWASLAB_EFFECT_ALLELE_COL).alias(COJO_A1_COL),
        nw.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(COJO_A2_COL),
        nw.col(GWASLAB_EFFECT_ALLELE_FREQ_COL).alias(COJO_FREQ_COL),
        nw.col(GWASLAB_BETA_COL).alias(COJO_B_COL),
        nw.col(GWASLAB_SE_COL).alias(COJO_SE_COL),
        nw.col(GWASLAB_P_COL).alias(COJO_P_COL),
    )
    with_n = _with_sample_size_column(renamed, phenotype)
    selected = with_n.select(COJO_COLUMN_ORDER)
    selected.collect().to_pandas().to_csv(out_path, index=False, sep="\t")
    return out_path
