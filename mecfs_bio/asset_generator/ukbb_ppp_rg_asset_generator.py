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
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.assets.reference_data.pqtls.raw.sun_et_al_2023_pqtls import (
    SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.ppp_ldsc.ppp_protein_rg_task import (
    PppProteinCrossTraitRgTask,
    PppRgConfig,
)
from mecfs_bio.constants.ppp_ldsc_constants import PPP_VARIANT_SET_CIS_EXCLUDED


@frozen
class PppRgTasks:
    rg_task: Task


def generate_ppp_rg_assets(
    base_name: str,
    trait_task: Task,
    config: PppRgConfig = PppRgConfig(),
    trait_pipe: DataProcessingPipe = IdentityPipe(),
):
    return PppRgTasks(
        rg_task=PppProteinCrossTraitRgTask.create(
            asset_id=base_name + "_rg_task",
            trait_task=trait_task,
            protein_tasks=HAPMAP_3_PPP_DATABASE.protein_tasks,
            index_task=HAPMAP_3_PPP_DATABASE_INDEX,
            ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
            sample_size_task=PPP_SAMPLE_SIZES,
            # Gene coordinates are only consulted to build cis windows.
            st3_task=(
                SUN_ET_AL_2023_PQTL_SUPPLEMENTARY_RAW
                if config.variant_set == PPP_VARIANT_SET_CIS_EXCLUDED
                else None
            ),
            config=config,
            trait_pipe=trait_pipe,
        ),
    )
