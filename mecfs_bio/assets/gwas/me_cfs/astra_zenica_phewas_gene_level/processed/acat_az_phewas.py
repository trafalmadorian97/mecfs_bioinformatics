from mecfs_bio.assets.gwas.me_cfs.astra_zenica_phewas_gene_level.raw.get_mecfs_az_phewas import (
    MECFS_AZ_PHEWAS,
)
from mecfs_bio.build_system.task.acat_task import AcatTask

MECFS_AZ_PHEWAS_ACAT = AcatTask.create(
    source_task=MECFS_AZ_PHEWAS,
    asset_id="mecfs_az_phewas_acat",
    group_by=["Gene", "Ancestry"],
    p_value_col="P value",
    model_col="Collapsing model",
    excluded_models=["syn"],
)
