"""
Genome-reference harmonization of Nalls et al. 2019's Parkinson's disease GWAS (without
23andMe; 15,056 cases + 18,618 UK Biobank proxy cases + 449,056 controls, European ancestry)
at the GWAS Catalog harmonised file's GRCh38 positions, against the UCSC hg38 genome and the
gnomAD v4.1 genomes allele-frequency panel (non-Finnish European frequencies).

The per-variant sample size is N_cases + N_controls. Indels are kept.
"""

import narwhals

from mecfs_bio.assets.gwas.parkinsons.nalls_et_al_2019.raw.nalls_2019_raw import (
    NALLS_2019_PARKINSONS_RAW,
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

NALLS_2019_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=NALLS_2019_PARKINSONS_RAW,
    target_asset_id=AssetId("nalls_2019_parkinsons_build_38_sumstats"),
    basic_check=True,
    genome_build="38",
    fmt=GWASLabColumnSpecifiers(
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        beta="beta",
        se="standard_error",
        p="p_value",
        eaf="effect_allele_frequency",
        n=GWASLAB_SAMPLE_SIZE_COLUMN,
    ),
    pre_pipe=CompositePipe(
        [
            DropNullsPipe(subset=["chromosome", "base_pair_location"]),
            ExprPipe(
                (narwhals.col("N_cases") + narwhals.col("N_controls")).alias(
                    GWASLAB_SAMPLE_SIZE_COLUMN
                )
            ),
        ]
    ),
)

NALLS_2019_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=NALLS_2019_BUILD_38_SUMSTATS,
        asset_id="nalls_2019_parkinsons_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

NALLS_2019_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="nalls_2019_parkinsons_build_38_genome_reference_harmonized",
    sumstats_task=NALLS_2019_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
    panel_ancestry="nfe",
)
