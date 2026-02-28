"""
The Syncope GWAS sumstats file is formatted in an usual way:
- Effect allele frequency is expressed as a percent, not a proportion.
- Some rows have multiple rsIDS, to account for rsid synonyms.

This Task preprocesses the file to address these issues, then creates GWASLAB sumstats
- Note 1: if the sumstats issues are not fixed prior to conversion to GWASlab sumstats format, this can cause downstream issues.
- Note 2:  This task resolves the issues of synonymous rsids by creating separate row for each synonym.  Downstream tasks tasks need to take this into account to avoid double-counting
"""

from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.raw.raw_syncope_data import (
    AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.scale_col_pipe import ScaleColPipe
from mecfs_bio.build_system.task.pipes.str_split_col import SplitColPipe
from mecfs_bio.build_system.task.pipes.unnest_pipe import UNNestPipe

SYNCOPE_SUMSTATS_EXPLODE_AND_SCALE = GWASLabCreateSumstatsTask(
    df_source_task=AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
    fmt=GWASLabColumnSpecifiers(
        chrom="Chrom",
        pos="Pos",
        ea="EA",
        nea="OA",
        eaf="EAF",
        p="Pval",
        snpid="ID",
        OR="OR",
        se="SE",
        info="Info",
        rsid="rsName",
    ),
    basic_check=True,
    genome_build="infer",
    pre_pipe=CompositePipe(
        [
            ScaleColPipe(
                col="EAF", scale_factor=1 / 100
            ),  # EAF is expressed as a percentage
            SplitColPipe(
                col_to_split="rsName", split_by=",", new_col_name="rsName"
            ),  # some rsid columns contain multiple synonymous rsids
            UNNestPipe(col_to_unnest="rsName"),
        ]
    ),
    asset_id=AssetId("aegisdottir_et_al_raw_sumstats_exploded_and_scaled"),
)
