import polars as pl
from pathlib import Path

import narwhals as nw

from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL

MIXER_RSID_COL  ="RSID"
MIXER_CHROM_COL = "CHR"
MIXER_POS_COL = "POS"
MIXER_EFFECT_ALLELE_COL = "EffectAllele"
MIXER_NON_EFFECT_ALLELE_COL = "OtherAllele"
MIXER_EFFECTIVE_SAMPLE_SIZE = "N"
MIXER_Z_SCORE_COL = "Z"

def _prep_summary_statistics_for_mixer(
    sumstats_dataframe_task:    Task,
    fetch: Fetch,
    temp_dir: Path,
    pipe: DataProcessingPipe
):
    asset = fetch(sumstats_dataframe_task.asset_id)
    frame= pipe.process(scan_dataframe_asset(asset, meta=sumstats_dataframe_task.meta)).collect().to_polars()
    frame.with_columns(
        pl.col(GWASLAB_RSID_COL).alias(MIXER_RSID_COL),
        pl.col(GWASLAB_CHROM_COL).alias(MIXER_CHROM_COL),
    pl.col(GWASLAB_POS_COL).alias(MIXER_POS_COL),
        pl.col
    )


