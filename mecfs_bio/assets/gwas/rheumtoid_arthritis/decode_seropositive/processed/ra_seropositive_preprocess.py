"""
Preprocessing step that
- Converts the whitespace-separated text-format summary statistics file to parquet, which is easier to work with
- Scale the effect allele frequency column so that it is expressed as a proportion, not a fraction.
- Restricts to variants with MAF>=0.05
"""

from mecfs_bio.assets.gwas.rheumtoid_arthritis.decode_seropositive.raw.download_raw_ra_seropositive import (
    DECODE_SEROPOSITIVE_RA_RAW,
)
from mecfs_bio.build_system.task.filter_snps_by_frequency import FilterSNPsFrequencyTask
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.scale_col_pipe import ScaleColPipe
from mecfs_bio.build_system.task.whitespace_sep_text_to_parquet_task import (
    WhitespaceSepTextToParquetTask,
)

PARQUET_SEROPOS_RA = WhitespaceSepTextToParquetTask.create(
    source_task=DECODE_SEROPOSITIVE_RA_RAW,
    asset_id="parquet_ra_seropositive",
)
SEROPOS_RA_EAF_SCALED = PipeDataFrameTask.create(
    source_task=PARQUET_SEROPOS_RA,
    asset_id="ra_seropos_eaf_scaled",
    out_format=ParquetOutFormat(),
    pipes=[ScaleColPipe(col="EAFrq", scale_factor=0.01)],
)
SEROPOS_RA_FILTERED_FOR_FREQ = FilterSNPsFrequencyTask.create(
    id=("ra_seropos_filter_for_freq"),
    raw_gwas_task=SEROPOS_RA_EAF_SCALED,
    allele_freq_col="EAFrq",
    freq_thresh=0.05,
)
