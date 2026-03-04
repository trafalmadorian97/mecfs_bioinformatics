from mecfs_bio.assets.gwas.ankylosing_spondylitis.ukbb.processed.ukbb_ank_spond_sumstats import \
    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GwasLabTransformSpec, HarmonizationOptions, \
    GWASLabVCFRef
from mecfs_bio.build_system.task.gwaslab.gwaslab_transform_sumstats import GWASLabTransformSumstatsTask

UKBB_ANK_SPOND_SUMSTATS_HARMONIZED = GWASLabTransformSumstatsTask.create_from_source_task(

    source_tsk=
    UK_BIOBANK_ANKYLOSING_SPONDYLITIS_SUMSTATS ,

    asset_id="uk_biobank_ank_spond_sumstats_37_harmonized",

    spec=GwasLabTransformSpec(
        harmonize_options=HarmonizationOptions(
            ref_infer=GWASLabVCFRef(name="1kg_eur_hg19", ref_alt_freq="AF"),
            ref_seq="ucsc_genome_hg19",
            check_ref_files=True,
            drop_missing_from_ref_seq=True,
            drop_missing_from_ref_infer_or_ambiguous=True,
            cores=1,
        )
    )

)
