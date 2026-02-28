"""
Task to restrict the Aegisdottir et al. summary statistics to hapmap3 and deduplicate them.
This avoids any potential issue downstream due to multiple rows referring to the same variant with different rsi
"""

from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.processed.syncope_sumstats_liftover import (
    AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GwasLabTransformSpec,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import (
    GWASLabTransformSumstatsTask,
)
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS_HAPMAP3_DEDUP = (
    GWASLabTransformSumstatsTask.create_from_source_task(
        source_tsk=AEGISDOTTIR_SYNCOPE_LIFTOVER_SUMSTATS,
        asset_id="aegissdotir_sumstats_liftover_hapmap3_dedup",
        spec=GwasLabTransformSpec(
            genome_build="infer", filter_hapmap3=True, liftover_to="19"
        ),
        post_pipe=UniquePipe(
            by=[
                GWASLAB_CHROM_COL,
                GWASLAB_POS_COL,
                GWASLAB_EFFECT_ALLELE_COL,
                GWASLAB_NON_EFFECT_ALLELE_COL,
            ],
            keep="first",
            order_by=[GWASLAB_RSID_COL],
        ),
    )
)
