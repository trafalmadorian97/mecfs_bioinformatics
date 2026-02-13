"""
Filter out columns referring to non-Finish Europeans


In the raw GWAS data, the NFE columns refer to non-Finish Europeans.
We extract these columns to get summary statistics for Non-Finish Europeans, and
drop any variants with null values in these columns, indicating they were not measured in the NFE population
"""

from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.raw_gwas_data.liu_et_al_2023_meta import (
    LIU_ET_AL_2023_IBD_META,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

LIU_ET_AL_2023_IBD_META_EUR_ONLY = PipeDataFrameTask.create(
    source_task=LIU_ET_AL_2023_IBD_META,
    asset_id="liu_et_al_2023_ibd_meta_eur_only",
    out_format=ParquetOutFormat(),
    pipes=[
        SelectColPipe(
            [
                "MarkerName",
                "Allele1",
                "Allele2",
                "CHR",
                "BP",
                "AF_NFE",
                "BETA_NFE",
                "SE_NFE",
                "P_NFE",
            ],
        ),
        DropNullsPipe(),
    ],
    backend="polars",
)
