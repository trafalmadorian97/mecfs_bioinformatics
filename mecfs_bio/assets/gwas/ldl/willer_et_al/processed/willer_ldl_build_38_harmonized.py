"""
Genome-reference harmonization of Willer et al.'s LDL GWAS at the GWAS Catalog's GRCh38
positions (hm_chrom, hm_pos), against the UCSC hg38 genome and the 1000 Genomes EUR (30x, GRCh38) allele-frequency panel.

The ORIGINAL alleles are used (effect_allele, other_allele, beta, standard_error), not the
Catalog's harmonised ones, and our own harmonizer orients them. The alleles are lowercase in the
file and are uppercased. eaf_ref is the 1000 Genomes EUR frequency of the original effect_allele,
except at palindromic SNVs (A/T, C/G), where it is a placeholder (null for about half of them
and exactly 1.0 for most of the rest); there it is set to null, so the untrusted harmonizer
drops them. Rows that failed the Catalog's liftover (no hm_pos) are dropped.
"""

import narwhals

from mecfs_bio.assets.gwas.ldl.willer_et_al.raw.raw_willer_ldl_data import (
    WILLER_LDL_EUR_DATA_RAW,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
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

_EA = "effect_allele"
_OA = "other_allele"
_EAF_COL = "eaf_ref_non_palindromic"
_UPPER_EA = narwhals.col(_EA).str.to_uppercase()
_UPPER_OA = narwhals.col(_OA).str.to_uppercase()
_IS_PALINDROMIC = (
    ((_UPPER_EA == "A") & (_UPPER_OA == "T"))
    | ((_UPPER_EA == "T") & (_UPPER_OA == "A"))
    | ((_UPPER_EA == "C") & (_UPPER_OA == "G"))
    | ((_UPPER_EA == "G") & (_UPPER_OA == "C"))
)

WILLER_LDL_BUILD_38_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=WILLER_LDL_EUR_DATA_RAW,
    target_asset_id=AssetId("willer_ldl_build_38_sumstats"),
    basic_check=True,
    genome_build="38",
    fmt=GWASLabColumnSpecifiers(
        chrom="hm_chrom",
        pos="hm_pos",
        ea=_EA,
        nea=_OA,
        beta="beta",
        se="standard_error",
        p="p_value",
        eaf=_EAF_COL,
        n="n",
    ),
    pre_pipe=CompositePipe(
        [
            DropNullsPipe(subset=["hm_chrom", "hm_pos"]),
            ExprPipe(
                narwhals.when(~_IS_PALINDROMIC)
                .then(narwhals.col("eaf_ref"))
                .alias(_EAF_COL)
            ),
            ExprPipe(_UPPER_EA.alias(_EA)),
            ExprPipe(_UPPER_OA.alias(_OA)),
        ]
    ),
)

WILLER_LDL_BUILD_38_PRE_HARMONIZATION_TABLE = (
    GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=WILLER_LDL_BUILD_38_SUMSTATS,
        asset_id="willer_ldl_build_38_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
    )
)

WILLER_LDL_BUILD_38_HARMONIZED = GenomeReferenceHarmonizationTask.create(
    asset_id="willer_ldl_build_38_genome_reference_harmonized",
    sumstats_task=WILLER_LDL_BUILD_38_PRE_HARMONIZATION_TABLE,
    fasta_task=UCSC_HG38_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
    panel_ancestry="eur",
)
