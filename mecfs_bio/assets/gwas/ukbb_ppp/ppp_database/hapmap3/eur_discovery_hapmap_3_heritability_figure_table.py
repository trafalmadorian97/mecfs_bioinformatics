"""
The per-protein HapMap3 EUR-discovery heritability table, encoded for display on
the documentation site by the data_table macro.

This is the same table as HAPMAP_3_PPP_HERITABILITY, re-encoded rather than
re-computed. At ~5,900 rows it is far too large for the markdown_table macro,
which renders every row into the page at build time; data_table instead ships
this parquet to the browser and renders it client-side.
"""

import narwhals

from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_heritability import (
    HAPMAP_3_PPP_HERITABILITY,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    ParquetWriteOptions,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    PipeDataFrameTask,
)
from mecfs_bio.build_system.task.pipes.attenuation_ratio_pipe import (
    AttenuationRatioPipe,
)
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.cast_to_float32_pipe import CastFloatsToFloat32Pipe
from mecfs_bio.build_system.task.pipes.compute_p_from_beta_se import (
    ComputePFromBetaSEPipeIfNeeded,
)
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.constants.ppp_ldsc_constants import (
    PPP_H2_H2_COL,
    PPP_H2_H2_SE_COL,
    PPP_H2_INTERCEPT_COL,
    PPP_H2_LAMBDA_GC_COL,
    PPP_H2_MEAN_CHI2_COL,
    PPP_H2_N_SNPS_COL,
)

HAPMAP_3_PPP_HERITABILITY_FIGURE_TABLE = PipeDataFrameTask.create(
    asset_id="ppp_heritability_hapmap_3_eur_discovery_table",
    source_task=HAPMAP_3_PPP_HERITABILITY,
    pipes=[
        CastPipe(
            target_column=PPP_H2_N_SNPS_COL,
            type=narwhals.Int32(),
            new_col_name=PPP_H2_N_SNPS_COL,
        ),
        ComputePFromBetaSEPipeIfNeeded(
            p_col="p",
            se_col=PPP_H2_H2_SE_COL,
            beta_col=PPP_H2_H2_COL,
        ),
        AttenuationRatioPipe(
            mean_chi_col=PPP_H2_MEAN_CHI2_COL, intercept_col=PPP_H2_INTERCEPT_COL
        ),
        DropColPipe([PPP_H2_LAMBDA_GC_COL]),
        CastFloatsToFloat32Pipe(),  # the input data for PPP LDSC is float32, so it is reasonable to cast the output to float32 as well
    ],
    out_format=ParquetOutFormat(
        write_options=ParquetWriteOptions(
            compression="zstd",
            compression_level=22,
            byte_stream_split_floats=True,
        )
    ),
    backend="polars",
)
