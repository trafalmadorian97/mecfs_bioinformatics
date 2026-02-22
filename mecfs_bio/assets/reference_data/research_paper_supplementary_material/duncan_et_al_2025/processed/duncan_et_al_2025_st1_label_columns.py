from mecfs_bio.assets.reference_data.research_paper_supplementary_material.duncan_et_al_2025.processed.duncan_et_al_2025_st1_extracted import (
    DUNCAN_ET_AL_2025_ST1_EXTRACTED,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    ParquetOutFormat,
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

DUNCAN_ET_AL_2025_ST1_LABEL_COLS = PipeDataFrameTask.create(
    source_task=DUNCAN_ET_AL_2025_ST1_EXTRACTED,
    asset_id="duncan_et_al_2025_st1_label_columns",
    out_format=ParquetOutFormat(),
    pipes=[
        SelectColPipe(
            [
                "VARIABLE",
                "Supercluster",
                "Class auto-annotation",
                "Neurotransmitter auto-annotation",
                "Neuropeptide auto-annotation",
                "Subtype auto-annotation",
                "Transferred MTG Label",
                "Top three regions",
                "Top Enriched Genes",
            ]
        )
    ],
)
