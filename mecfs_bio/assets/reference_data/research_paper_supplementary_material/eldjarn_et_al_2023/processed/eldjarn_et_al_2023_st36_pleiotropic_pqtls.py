"""
The highly pleiotropic pQTLs of Eldjarn et al. 2023 ST36 with typed, renamed columns (see
St36PleiotropicPqtlsPipe), for flagging deCODE coloc results near them.
"""

from mecfs_bio.assets.reference_data.research_paper_supplementary_material.eldjarn_et_al_2023.processed.eldjarn_et_al_2023_st36_extracted import (
    ELDJARN_ET_AL_2023_ST36_EXTRACTED,
)
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.st36_pleiotropic_pqtls_pipe import (
    St36PleiotropicPqtlsPipe,
)

ELDJARN_ET_AL_2023_ST36_PLEIOTROPIC_PQTLS = PipeDataFrameTask.create(
    source_task=ELDJARN_ET_AL_2023_ST36_EXTRACTED,
    asset_id="eldjarn_et_al_2023_st36_pleiotropic_pqtls",
    out_format=ParquetOutFormat(),
    pipes=[St36PleiotropicPqtlsPipe()],
    backend="polars",
)
