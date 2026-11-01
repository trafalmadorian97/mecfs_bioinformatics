"""
Genome-reference harmonization of deCODE's seropositive RA meta-GWAS at its GRCh38 positions
(PosB38), against the UCSC hg38 genome and the 1000 Genomes EUR (30x, GRCh38) allele-frequency panel.

Starts from the EAF-rescaled table restricted to MAF >= 0.01: the full table has 61M variants,
which gwaslab cannot hold on a 16 GB machine, while MAF >= 0.01 keeps 14.2M, a density
comparable to the UKB-PPP imputed set. (The standard analysis's MAF >= 0.05 filter would drop
variants colocalization wants.) As in the standard analysis: only SNVs are kept, the per-variant effective N comes
from the Cohorts string, and BETA and SE are computed from OR and P (after gwaslab, on its
column names).
"""

from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.processed.ra_seropositive_preprocess import (
    SEROPOS_RA_EAF_SCALED,
)
from mecfs_bio.assets.gwas.rheumtoid_arthritis.util.decode_cohort_sample_sizes import (
    load_decode_ra_cohorts,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_beta_pipe import ComputeBetaPipe
from mecfs_bio.build_system.task.pipes.compute_se_pipe import ComputeSEPipe
from mecfs_bio.build_system.task.pipes.effective_n_from_cohort_string_pipe import (
    EffectiveNFromCohortStringPipe,
)
from mecfs_bio.build_system.task.pipes.filter_rows_by_value import FilterRowsByValue
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

_ACGT = ["A", "C", "G", "T"]

SEROPOS_RA_FILTERED_FOR_FREQ_0_01 = FilterSNPsFrequencyTask.create(
    id="ra_seropos_filter_for_freq_0_01",
    raw_gwas_task=SEROPOS_RA_EAF_SCALED,
    allele_freq_col="EAFrq",
    freq_thresh=0.01,
)

SEROPOS_RA_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=SEROPOS_RA_FILTERED_FOR_FREQ_0_01,
    target_asset_id=AssetId("ra_seropositive_build_38_sumstats"),
    basic_check=True,
    genome_build="38",
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
    pre_pipe=CompositePipe(
        [
            FilterRowsByValue(target_column="EA", valid_values=_ACGT),
            FilterRowsByValue(target_column="OA", valid_values=_ACGT),
            EffectiveNFromCohortStringPipe(
                cohorts=load_decode_ra_cohorts("seropositive"),
                cohort_string_col="Cohorts",
                out_col=GWASLAB_SAMPLE_SIZE_COLUMN,
            ),
        ]
    ),
)

SEROPOS_RA_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=SEROPOS_RA_BUILD_38_SUMSTATS,
        asset_id="ra_seropositive_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
        pipe=CompositePipe([ComputeBetaPipe(), ComputeSEPipe()]),
    )
)

SEROPOS_RA_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="ra_seropositive_build_38_genome_reference_harmonized",
    sumstats_task=SEROPOS_RA_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
    panel_ancestry="eur",
)
