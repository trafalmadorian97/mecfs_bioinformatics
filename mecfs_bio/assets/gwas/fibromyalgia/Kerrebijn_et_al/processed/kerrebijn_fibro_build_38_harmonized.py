"""
Genome-reference harmonization of Kerrebijn et al.'s European fibromyalgia GWAS at the GWAS
Catalog harmonised file's GRCh38 positions, against the UCSC hg38 genome and the gnomAD v4.1
genomes allele-frequency panel (non-Finnish European frequencies).

The per-variant sample size is the deposited n_eff (the effective number of a balanced
case-control study), so the coloc trait model uses case fraction 0.5.
"""

import narwhals

from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.raw.kerrebijin_et_al_fibro_raw import (
    KERREBIJN_ET_AL_FIBRO_EUR_RAW,
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
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe

KERREBIJN_FIBRO_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=KERREBIJN_ET_AL_FIBRO_EUR_RAW,
    target_asset_id=AssetId("kerrebijn_fibro_build_38_sumstats"),
    basic_check=True,
    genome_build="38",
    fmt=GWASLabColumnSpecifiers(
        rsid="rsid",
        chrom="chromosome",
        pos="base_pair_location",
        ea="effect_allele",
        nea="other_allele",
        beta="beta",
        se="standard_error",
        p="p_value",
        eaf="effect_allele_frequency",
        n="n_eff",
    ),
    pre_pipe=CastPipe(
        target_column="n_eff", type=narwhals.Int64(), new_col_name="n_eff"
    ),
)

KERREBIJN_FIBRO_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=KERREBIJN_FIBRO_BUILD_38_SUMSTATS,
        asset_id="kerrebijn_fibro_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

KERREBIJN_FIBRO_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="kerrebijn_fibro_build_38_genome_reference_harmonized",
    sumstats_task=KERREBIJN_FIBRO_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
    panel_ancestry="nfe",
)
