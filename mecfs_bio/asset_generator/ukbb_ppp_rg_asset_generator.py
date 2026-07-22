from attrs import frozen

from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_database_protein_files import (
    HAPMAP_3_PPP_DATABASE,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.eur_discovery_hapmap3_ppp_heritability import (
    PPP_SAMPLE_SIZES,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genomes_phase_3_v1_consolidated import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
)
from mecfs_bio.assets.reference_data.pqtls.processed.sun_et_al_2023_st3_extracted import (
    SUN_ET_AL_2023_ST3_EXTRACTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat, ParquetWriteOptions
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.cast_to_float32_pipe import CastFloatsToFloat32Pipe
from mecfs_bio.build_system.task.pipes.cast_to_int32pipe import CastIntsToInt32Pipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.drop_nan_pipe import DropNanPipe
from mecfs_bio.build_system.task.pipes.filter_rows_by_min_in_col import FilterRowsByMinInCol
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import (
    PppProteinCrossTraitRgTask,
    PppRgConfig,
)
from mecfs_bio.constants.ppp_ldsc_constants import PPP_VARIANT_SET_CIS_EXCLUDED, PPP_RG_RG_Z_COL, PPP_RG_H2_TRAIT_COL, \
    PPP_RG_N_BAR_PROTEIN_COL, PPP_RG_N_BAR_TRAIT_COL, PPP_RG_RG_COL, PPP_RG_RG_SE_COL, PPP_RG_H2_PROTEIN_COL


@frozen
class PppRgTasks:
    rg_task: Task
    display_frame_task: Task

    def get_terminal_tasks(self) -> list[Task]:
        return [self.rg_task, self.display_frame_task]


def generate_ppp_rg_assets(
    base_name: str,
    trait_task: Task,
    config: PppRgConfig = PppRgConfig(),
    trait_pipe: DataProcessingPipe = IdentityPipe(),
):
    rg_task = PppProteinCrossTraitRgTask.create(
            asset_id=base_name + "_rg_task",
            trait_task=trait_task,
            protein_tasks=HAPMAP_3_PPP_DATABASE.protein_tasks,
            index_task=HAPMAP_3_PPP_DATABASE_INDEX,
            ld_scores_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_CONSOLIDATE,
            sample_size_task=PPP_SAMPLE_SIZES,
            # Gene coordinates are only consulted to build cis windows.
            gene_coords_task=(
                SUN_ET_AL_2023_ST3_EXTRACTED
                if config.variant_set == PPP_VARIANT_SET_CIS_EXCLUDED
                else None
            ),
            config=config,
            trait_pipe=trait_pipe,
        )

    display_frame_task = PipeDataFrameTask.create(
        source_task=rg_task,
        pipes=[
            DropColPipe(
                [
                    PPP_RG_RG_Z_COL,
                    PPP_RG_N_BAR_PROTEIN_COL,
                    PPP_RG_N_BAR_TRAIT_COL


                ]
            ),
            DropNanPipe(cols=[
                PPP_RG_RG_COL,
                PPP_RG_RG_SE_COL
            ]),
            FilterRowsByMinInCol(
              min_value=0.02,
                col=PPP_RG_H2_PROTEIN_COL
            ),
            CastFloatsToFloat32Pipe(),
            CastIntsToInt32Pipe(),
        ],
        asset_id=base_name + "_display_frame",
        out_format=ParquetOutFormat(
            ParquetWriteOptions(
                compression="zstd",
                compression_level=22,
                byte_stream_split_floats=True
            )
        )

    )

    return PppRgTasks(
        rg_task=rg_task,
        display_frame_task=display_frame_task,
    )
