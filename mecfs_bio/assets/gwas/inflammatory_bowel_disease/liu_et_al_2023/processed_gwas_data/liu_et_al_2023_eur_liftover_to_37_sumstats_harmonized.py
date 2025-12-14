from mecfs_bio.assets.gwas.inflammatory_bowel_disease.liu_et_al_2023.processed_gwas_data.liu_et_al_2023_eur_liftover_to_37_sumstats import (
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GwasLabTransformSpec,
    GWASLabVCFRef,
    HarmonizationOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import (
    GWASLabTransformSumstatsTask,
)

LIU_ET_AL_2023_IBD_EUR_HARMONIZE = GWASLabTransformSumstatsTask.create_from_source_task(
    LIU_ET_AL_2023_IBD_EUR_LIFTOVER_37_SUMSTATS,
    asset_id="liu_et_al_2023_ibd_37_harmonized",
    spec=GwasLabTransformSpec(
        harmonize_options=HarmonizationOptions(
            ref_infer=GWASLabVCFRef(name="1kg_eur_hg19", ref_alt_freq="AF"),
            ref_seq="ucsc_genome_hg19",
            check_ref_files=True,
            drop_missing_from_ref=True,
            cores=4,
        )
    ),
)
