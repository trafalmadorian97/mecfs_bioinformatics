from mecfs_bio.asset_generator.h_magma_asset_generator import generate_h_magma_tasks
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import \
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

DECODE_ME_H_MAGMA_ASSET_GENERATOR =generate_h_magma_tasks(
    base_name="decode_me",
    gwas_parquet_with_rsids_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
    sample_size=275488, # from preprint
    pipes=[
RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL)
    ],
)