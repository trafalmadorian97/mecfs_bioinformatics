from mecfs_bio.assets.gwas.me_cfs.astra_zenica_phewas_gene_level.raw.get_mecfs_az_phewas import (
    MECFS_AZ_PHEWAS,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.pivot_pipe import PivotPipe

MECFS_AZ_PHEWAS_PIVOTED_ON_P = PipeDataFrameTask.create(
    MECFS_AZ_PHEWAS,
    asset_id="mecfs_az_phewas_pivoted_on_p",
    out_format=ParquetOutFormat(),
    backend="polars",
    pipes=[
        PivotPipe(
            index=["Gene", "Ancestry"],
            on="Collapsing model",
            values="P value",
        )
    ],
)
