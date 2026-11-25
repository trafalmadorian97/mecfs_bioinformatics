"""
Genome-reference harmonization of the PGC 2022 schizophrenia GWAS (Trubetskoy et al., European
autosomal summary statistics, 53,386 cases and 77,258 controls), lifted over from GRCh37 to
GRCh38 by gwaslab and harmonized against the UCSC hg38 genome and the gnomAD v4.1 genomes
allele-frequency panel (non-Finnish European frequencies).

A1 is the effect allele and BETA its log odds ratio. The effect allele frequency is FCON, the
frequency in controls (closest to a population frequency). The per-variant sample size is
NCAS + NCON. Variants gwaslab cannot lift over are dropped. Indels are kept.
"""

import narwhals

from mecfs_bio.assets.gwas.schizophrenia.pgc2022.raw.raw_sch_pgc2022 import (
    PGC_2022_SCH_RAW,
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

PGC_2022_SCH_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=PGC_2022_SCH_RAW,
    target_asset_id=AssetId("pgc_2022_sch_build_38_sumstats"),
    basic_check=True,
    genome_build="19",
    liftover_to="38",
    fmt=GWASLabColumnSpecifiers(
        rsid="ID",
        chrom="CHROM",
        pos="POS",
        ea="A1",
        nea="A2",
        beta="BETA",
        se="SE",
        p="PVAL",
        eaf="FCON",
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
    pre_pipe=ExprPipe(
        (narwhals.col("NCAS") + narwhals.col("NCON")).alias(GWASLAB_SAMPLE_SIZE_COLUMN)
    ),
)

PGC_2022_SCH_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=PGC_2022_SCH_BUILD_38_SUMSTATS,
        asset_id="pgc_2022_sch_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

PGC_2022_SCH_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="pgc_2022_sch_build_38_genome_reference_harmonized",
    sumstats_task=PGC_2022_SCH_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
    panel_ancestry="nfe",
)
