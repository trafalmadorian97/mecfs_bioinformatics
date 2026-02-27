from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_explode_and_scale import \
    SYNCOPE_SUMSTATS_EXPLODE_AND_SCALE
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GwasLabTransformSpec
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import GWASLabTransformSumstatsTask

AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS =GWASLabTransformSumstatsTask.create_from_source_task(
    source_tsk=SYNCOPE_SUMSTATS_EXPLODE_AND_SCALE,
    asset_id="aegissdotir_sumstats_liftover",
    spec=GwasLabTransformSpec(
        genome_build="infer",
        liftover_to="19"
    )
)