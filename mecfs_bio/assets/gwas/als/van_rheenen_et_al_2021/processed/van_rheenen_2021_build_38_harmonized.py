"""
Genome-reference harmonization of van Rheenen et al. 2021's European-ancestry ALS GWAS, lifted
over from GRCh37 to GRCh38 by gwaslab and harmonized against the UCSC hg38 genome and the
gnomAD v4.1 genomes allele-frequency panel (non-Finnish European frequencies).

Alleles are uppercased. The per-variant sample size is the deposited N_effective (the
effective number of a balanced case-control study), so the coloc trait model uses case
fraction 0.5. Variants gwaslab cannot lift over are dropped. Indels are kept.
"""

import narwhals

from mecfs_bio.assets.gwas.als.van_rheenen_et_al_2021.raw.van_rheenen_2021_raw import (
    VAN_RHEENEN_2021_ALS_EUR_RAW,
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
from mecfs_bio.build_system.task.pipes.expr_pipe import ExprPipe

_EA = "effect_allele"
_OA = "other_allele"

VAN_RHEENEN_2021_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=VAN_RHEENEN_2021_ALS_EUR_RAW,
    target_asset_id=AssetId("van_rheenen_2021_als_eur_build_38_sumstats"),
    basic_check=True,
    genome_build="19",
    liftover_to="38",
    fmt=GWASLabColumnSpecifiers(
        rsid="rsid",
        chrom="chromosome",
        pos="base_pair_location",
        ea=_EA,
        nea=_OA,
        beta="beta",
        se="standard_error",
        p="p_value",
        eaf="effect_allele_frequency",
        n="N_effective",
    ),
    pre_pipe=CompositePipe(
        [
            ExprPipe(narwhals.col(_EA).str.to_uppercase().alias(_EA)),
            ExprPipe(narwhals.col(_OA).str.to_uppercase().alias(_OA)),
        ]
    ),
)

VAN_RHEENEN_2021_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=VAN_RHEENEN_2021_BUILD_38_SUMSTATS,
        asset_id="van_rheenen_2021_als_eur_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

VAN_RHEENEN_2021_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="van_rheenen_2021_als_eur_build_38_genome_reference_harmonized",
    sumstats_task=VAN_RHEENEN_2021_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=GNOMAD_V4_1_GENOMES_HG38_ALLELE_FREQUENCIES,
    panel_ancestry="nfe",
)
