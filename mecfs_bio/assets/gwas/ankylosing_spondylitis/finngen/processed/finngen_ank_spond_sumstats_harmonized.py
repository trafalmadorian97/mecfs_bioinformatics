from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.processed.finngen_ank_spond_sumstats import (
    FINGNEN_ANK_SPOND_SUMSTATS,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GwasLabTransformSpec,
    GWASLabVCFRef,
    HarmonizationOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import (
    GWASLabTransformSumstatsTask,
)

FINGNEN_ANK_SPOND_SUMSTATS_HARMONIZED = (
    GWASLabTransformSumstatsTask.create_from_source_task(
        source_tsk=FINGNEN_ANK_SPOND_SUMSTATS,
        asset_id="finngen_ank_spond_sumstats_37_harmonized",
        spec=GwasLabTransformSpec(
            harmonize_options=HarmonizationOptions(
                ref_infer=GWASLabVCFRef(name="1kg_eur_hg19", ref_alt_freq="AF"),
                ref_seq="ucsc_genome_hg19",
                check_ref_files=True,
                drop_missing_from_ref_seq=True,
                drop_missing_from_ref_infer_or_ambiguous=True,
                cores=1,
            )
        ),
    )
)
