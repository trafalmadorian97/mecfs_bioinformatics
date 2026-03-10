from mecfs_bio.assets.reference_data.lava_locus_file.default.raw.default_lava_locus_file import DEFAULT_LAVA_LOCUS_FILE
from mecfs_bio.build_system.task.fuse_loci_for_lava_task import FuseLociForLavaTask

DEFAULT_LAVA_LOCUS_FUSED = FuseLociForLavaTask.create(
    asset_id="default_lava_locus_fused",
    original_loci_task=DEFAULT_LAVA_LOCUS_FILE,
    original_loci_per_fused_locus=10
)