import narwhals

from mecfs_bio.asset_generator.concrete_standard_analysis_task_generator import (
    concrete_standard_analysis_generator_assume_already_has_rsid,
)
from mecfs_bio.assets.gwas.fibromyalgia.Kerrebijn_et_al.raw.kerrebijin_et_al_fibro_raw import (
    KERREBIJN_ET_AL_FIBRO_EUR_RAW,
)
from mecfs_bio.build_system.sample_size_spec import PerVariantSampleSize
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe

KERREBIJN_ET_AL_FIBRO_STANDARD_ANALYSIS = (
    concrete_standard_analysis_generator_assume_already_has_rsid(
        base_name="kerrebijin_fibro",
        raw_gwas_data_task=KERREBIJN_ET_AL_FIBRO_EUR_RAW,
        fmt=GWASLabColumnSpecifiers(
            chrom="chromosome",
            pos="base_pair_location",
            ea="effect_allele",
            nea="other_allele",
            beta="beta",
            se="standard_error",
            eaf="effect_allele_frequency",
            p="p_value",
            rsid="rsid",
            n="n_eff",  # use effective n so we can set case/control to 1:1
        ),
        sample_size=PerVariantSampleSize(),
        pre_pipe=CastPipe(
            target_column="n_eff", type=narwhals.Int64(), new_col_name="n_eff"
        ),
    )
)
