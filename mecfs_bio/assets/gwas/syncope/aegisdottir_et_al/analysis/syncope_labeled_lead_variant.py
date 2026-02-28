from mecfs_bio.asset_generator.labeled_lead_variants_asset_generator import (
    generate_tasks_labeled_lead_variants,
)
from mecfs_bio.assets.gwas.syncope.aegisdottir_et_al.raw.raw_syncope_data import (
    AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabColumnSpecifiers,
)

SYNCOPE_LABELED_LEAD_VARIANT = generate_tasks_labeled_lead_variants(
    base_name="syncope_labeled_lead_variant_tasks",
    raw_gwas_data_task=AEGISDOTTIR_ET_AL_RAW_SYNCOPE_DATA,
    fmt=GWASLabColumnSpecifiers(
        chrom="Chrom",
        pos="Pos",
        ea="EA",
        nea="OA",
        eaf="EAF",
        p="Pval",
        snpid="ID",
        OR="OR",
        se="SE",
        info="Info",
        rsid="rsName",
    ),
)
