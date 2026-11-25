"""
Genome-reference harmonization of Bellenguez et al. 2022's Alzheimer's disease GWAS at the GWAS
Catalog's GRCh38 harmonised positions (hm_chrom, hm_pos), against the UCSC hg38 genome and the
gnomAD v4.1 genomes allele-frequency panel (non-Finnish European frequencies).

The Catalog's harmonised alleles, beta and frequency are used (hm_effect_allele,
hm_other_allele, hm_beta, hm_effect_allele_frequency); our own harmonizer orients them again.
The per-variant sample size is n_cas + n_con, cases including the UK Biobank proxy cases
(family history of dementia). Rows without a harmonised position are dropped. Indels are kept.
"""

import narwhals

from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.raw.raw_bellenguez_data import (
    BELLENGUEZ_ET_AL_ALZHEIMERS_RAW,
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
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.expr_pipe import ExprPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_SAMPLE_SIZE_COLUMN

BELLENGUEZ_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=BELLENGUEZ_ET_AL_ALZHEIMERS_RAW,
    target_asset_id=AssetId("bellenguez_alzheimers_build_38_sumstats"),
    basic_check=True,
    genome_build="38",
    fmt=GWASLabColumnSpecifiers(
        chrom="hm_chrom",
        pos="hm_pos",
        ea="hm_effect_allele",
        nea="hm_other_allele",
        beta="hm_beta",
        se="standard_error",
        p="p_value",
        eaf="hm_effect_allele_frequency",
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
    pre_pipe=CompositePipe(
        [
            DropNullsPipe(subset=["hm_chrom", "hm_pos"]),
            ExprPipe(
                (narwhals.col("n_cas") + narwhals.col("n_con")).alias(
                    GWASLAB_SAMPLE_SIZE_COLUMN
                )
            ),
        ]
    ),
)

BELLENGUEZ_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=BELLENGUEZ_BUILD_38_SUMSTATS,
        asset_id="bellenguez_alzheimers_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

BELLENGUEZ_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="bellenguez_alzheimers_build_38_genome_reference_harmonized",
    sumstats_task=BELLENGUEZ_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
    panel_ancestry="nfe",
)
