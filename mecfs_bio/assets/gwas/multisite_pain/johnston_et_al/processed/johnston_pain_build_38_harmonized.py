"""
Genome-reference harmonization of Johnston et al.'s UK Biobank multisite chronic pain GWAS
(BOLT-LMM, quantitative: number of chronic pain sites), lifted over from GRCh37 to GRCh38 by
gwaslab and harmonized against the UCSC hg38 genome and the gnomAD v4.1 genomes
allele-frequency panel (non-Finnish European frequencies).

ALLELE1 is the effect allele, with BETA and SE from the BOLT-LMM infinitesimal model and
P_BOLT_LMM_INF as the p value. The file has no sample-size column, so every variant gets the
GWAS Catalog total of 387,649. Variants gwaslab cannot lift over are dropped.
"""

import narwhals

from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.auxiliary.sample_info import (
    JOHNSTON_ET_AL_SAMPLE_INFO,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.processed.extracted_johnston_data import (
    JOHNSTON_ET_AL_PAIN_EXTRACTED,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
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
from mecfs_bio.build_system.task.pipes.expr_pipe import ExprPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

_TOTAL_SAMPLE_SIZE = JOHNSTON_ET_AL_SAMPLE_INFO.total_sample_size
assert _TOTAL_SAMPLE_SIZE is not None

JOHNSTON_PAIN_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=JOHNSTON_ET_AL_PAIN_EXTRACTED,
    target_asset_id=AssetId("johnston_et_al_pain_build_38_sumstats"),
    basic_check=True,
    genome_build="19",
    liftover_to="38",
    fmt=GWASLabColumnSpecifiers(
        rsid="SNP",
        chrom="CHR",
        pos="BP",
        ea="ALLELE1",
        nea="ALLELE0",
        beta="BETA",
        se="SE",
        p="P_BOLT_LMM_INF",
        eaf="A1FREQ",
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
    pre_pipe=ExprPipe(
        narwhals.lit(_TOTAL_SAMPLE_SIZE).alias(GWASLAB_SAMPLE_SIZE_COLUMN)
    ),
)

JOHNSTON_PAIN_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=JOHNSTON_PAIN_BUILD_38_SUMSTATS,
        asset_id="johnston_et_al_pain_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

JOHNSTON_PAIN_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="johnston_et_al_pain_build_38_genome_reference_harmonized",
    sumstats_task=JOHNSTON_PAIN_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
    panel_ancestry="nfe",
)
