from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_no_rsid,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.raw.download_raw_ra_seropositive import (
    DECODE_SEROPOSITIVE_RA_RAW,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.util.decode_cohort_sample_sizes import (
    load_decode_ra_cohorts,
)
from mecfs_bio.build_system.sample_size_spec import PerVariantSampleSize
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.effective_n_from_cohort_string_pipe import (
    EffectiveNFromCohortStringPipe,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.build_system.task.whitespace_sep_text_to_parquet_task import (
    WhitespaceSepTextToParquetTask,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

# The deCODE summary statistics have no sample-size column; instead each variant's
# Cohorts string reports which of the six contributing cohorts backed it. Derive a
# per-variant effective sample size from that string before sumstats creation, so
# MAGMA, S-LDSC, LDSC heritability and ct-LDSC genetic correlation all read N from
# the data rather than a single constant.
_COHORT_STRING_COL = "Cohorts"

_add_effective_n = EffectiveNFromCohortStringPipe(
    cohorts=load_decode_ra_cohorts("seropositive"),
    cohort_string_col=_COHORT_STRING_COL,
    out_col=GWASLAB_SAMPLE_SIZE_COLUMN,
)

# Restrict to biallelic SNPs before gwaslab. deCODE retains multi-base alleles (MNVs
# and long substitutions up to ~800bp) that survive gwaslab's indel filter, because
# that filter only drops length-mismatched indels. Harmonization's reverse-complement
# then runs Categorical.astype(str) on the EA/NEA columns, allocating
# rows x max-allele-length x 4B (~50GB at 16M rows x 814 chars) -> OOM. Keeping only
# single-base A/C/T/G alleles makes every allele length 1, and is the variant set that
# MAGMA / S-LDSC / ct-LDSC use anyway. It also lets gwaslab skip the slow indel
# normalization step.
_ACTG = ["A", "C", "T", "G"]
_pre_gwaslab_pipe = CompositePipe(
    [
        FilterRowsByValue(target_column="EA", valid_values=_ACTG),
        FilterRowsByValue(target_column="OA", valid_values=_ACTG),
        _add_effective_n,
    ]
)

# The deCODE file is arbitrary-whitespace-separated, which polars cannot scan, so the
# read_spec machinery would fall back to a single eager pandas read of all 61M rows
# (OOM). Convert to parquet once in bounded-memory chunks; every downstream task then
# scans it lazily with projection pushdown.
PARQUET_SEROPOS_RA = WhitespaceSepTextToParquetTask.create(
    source_task=DECODE_SEROPOSITIVE_RA_RAW,
    asset_id="parquet_ra_seropositive",
)

SCALED_PARQUET

SEROPOS_RA_FILTERED_FOR_FREQ = (
    FilterSNPsFrequencyTask.create(
    id=("ra_seropos_filter_for_freq"),
    raw_gwas_task=PARQUET_SEROPOS_RA,
    allele_freq_col="EAFrq",
    freq_thresh=0.05,
))

DECODE_SEROPOSITIVE_RA_STANDARD_ANALYSIS = concrete_standard_analysis_generator_no_rsid(
    base_name="decode_ra_seropositive",
    raw_gwas_data_task=SEROPOS_RA_FILTERED_FOR_FREQ,  # PARQUET_SEROPOS_RA,
    fmt=GWASLabColumnSpecifiers(
        chrom="Chr",
        pos="PosB38",
        # Deliberately omit snpid: the deCODE "chr:pos_ref_alt" ID is the single
        # largest string column and would dominate the pandas materialisation in
        # gwaslab sumstats creation (~24GB peak at 45.8M variants -> OOM on a 16GB
        # box). rsids are reassigned downstream via annovar, so gwaslab can
        # regenerate SNPID from chrom:pos:ea:nea; dropping it halves the peak.
        ea="EA",
        nea="OA",
        eaf="EAFrq",
        OR="OR",
        p="P",
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
    sample_size=PerVariantSampleSize(),
    pre_pipe_before_rsid_assignment=_pre_gwaslab_pipe,
    # deCODE reports OR and P but no BETA/SE. Compute BETA = ln(OR), then SE from
    # BETA and P, so S-LDSC / cell-type-specific LDSC can form Z = BETA / SE (it
    # errors otherwise: "Cannot create Z: need either Z column, or BETA and SE").
    pre_pipe_after_rsid_assignment=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    # Redundant safety net now that _pre_gwaslab_pipe keeps only SNPs, but harmless and
    # documents that no indels should reach harmonization.
    filter_indels_in_harmonized=True,
)
