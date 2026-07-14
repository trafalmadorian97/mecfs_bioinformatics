"""
Task generator to run standard analysis on DECODE's meta-gwas of rheumatoid arthritis.
Additional preprocessing steps:
- Filter out non-SNP variants for simplicity, and to save memory. GWASLab 4.1.9 has a memory issue when running harmonization on long indels
- Compute effective N for each variant, based on the column that indicates which cohorts are used to compute which variants.
"""

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_no_rsid,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seronegative.processed.ra_seronegative_preprocess import (
    SERONEG_RA_FILTERED_FOR_FREQ,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.util.decode_cohort_sample_sizes import (
    load_decode_ra_cohorts,
)
from mecfs_bio.build_system.sample_size_spec import PerVariantSampleSize
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.effective_n_from_cohort_string_pipe import (
    EffectiveNFromCohortStringPipe,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

_COHORT_STRING_COL = "Cohorts"
_add_effective_n = EffectiveNFromCohortStringPipe(
    cohorts=load_decode_ra_cohorts("seronegative"),
    cohort_string_col=_COHORT_STRING_COL,
    out_col=GWASLAB_SAMPLE_SIZE_COLUMN,
)
_ACTG = ["A", "C", "T", "G"]
_pre_gwaslab_pipe = CompositePipe(
    [
        FilterRowsByValue(target_column="EA", valid_values=_ACTG),
        FilterRowsByValue(target_column="OA", valid_values=_ACTG),
        _add_effective_n,
    ]
)
SERONEGATIVE_RA_STANDARD_ANALYSIS = concrete_standard_analysis_generator_no_rsid(
    base_name="decode_ra_seronegative",
    raw_gwas_data_task=SERONEG_RA_FILTERED_FOR_FREQ,  # PARQUET_SERONEG_RA,
    fmt=GWASLabColumnSpecifiers(
        chrom="Chr",
        pos="PosB38",
        ea="EA",
        nea="OA",
        eaf="EAFrq",
        OR="OR",
        p="P",
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
    sample_size=PerVariantSampleSize(),
    pre_pipe_before_rsid_assignment=_pre_gwaslab_pipe,
    pre_pipe_after_rsid_assignment=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    filter_indels_in_harmonized=True,
    phenotype_info_for_ldsc=BinaryPhenotypeSampleInfo(
        sample_prevalence=0.5, estimated_population_prevalence=0.005
    ),  # Sample prevalence must be set to 0.5 since we are working with effective sample sizes.
    # Pop prevalence comes from https://pmc.ncbi.nlm.nih.gov/articles/PMC8053349/ which estimates 1% for RA in general.  Divide by 2 for a rough estimate for seronegative RA.
)
