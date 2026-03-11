"""
Initial exploratory LAVA run.  Will add more comprehensive LAVA runs later.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.auxiliary.prevalance_info import (
    DECODE_ME_PREVALENCE_INFO,
)
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS,
)
from mecfs_bio.assets.gwas.multi_trait.genetic_correlation.ct_ldsc.ct_ldsc_initial_asset_generator import (
    CT_LDSC_INITIAL_ASSET_GENERATOR,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.lava.lava_ld_reference.g1000_eur.processed.lava_thousand_geomes_eur_ld_ref_extracted import (
    LAVA_G100_EUR_LD_REF_EXTRACTED,
)
from mecfs_bio.assets.reference_data.lava.lava_locus_file.default.raw.default_lava_locus_file import (
    DEFAULT_LAVA_LOCUS_FILE,
)
from mecfs_bio.build_system.task.lava_task import (
    LavaBinarySampleInfo,
    LavaPhenotypeDataSource,
    LavaTask,
    LDReferenceInfo,
)
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.set_col_pipe import SetColToConstantPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
)

BASIC_G100_LAVA_ANALYSIS = LavaTask.create(
    asset_id="initial_g1000_lava_analysis",
    sources=[
        LavaPhenotypeDataSource(
            task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
            alias="DecodeME",
            sample_info=LavaBinarySampleInfo.from_ct_ldsc_sample_info(
                DECODE_ME_PREVALENCE_INFO
            ),
            pipe=RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
        ),
        LavaPhenotypeDataSource(
            JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            alias="Multisite_pain",
            pipe=SetColToConstantPipe(GWASLAB_SAMPLE_SIZE_COLUMN, constant=387649),
        ),
    ],
    ld_reference_info=LDReferenceInfo(
        ld_ref_task=LAVA_G100_EUR_LD_REF_EXTRACTED,
        filename_prefix="g1000_eur",
    ),
    lava_locus_definitions_task=DEFAULT_LAVA_LOCUS_FILE,
    ct_ldsc_task_for_overlap=CT_LDSC_INITIAL_ASSET_GENERATOR.aggregation_task,
    heritability_task_for_overlap=CT_LDSC_INITIAL_ASSET_GENERATOR.heritability_aggregation_task,
    max_loci=10,
)
