from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.raw.eur_ld_scores_thousand_genome_phase_3_v1_raw import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_RAW,
)
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask

THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED = (
    ExtractTarGzipTask.create(
        asset_id="thousand_genomes_phase_3_v1_eur_ld_scores_extracted",
        source_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_RAW,
        sub_folder_name_inside_tar="LDscore",
    )
)
