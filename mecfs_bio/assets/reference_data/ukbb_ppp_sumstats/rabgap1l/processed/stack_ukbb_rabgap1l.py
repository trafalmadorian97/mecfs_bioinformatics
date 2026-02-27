from mecfs_bio.assets.reference_data.ukbb_ppp_sumstats.rabgap1l.processed.untar_ukbb_rabgap_1l import (
    UKBB_PPP_RABGAP1L_UNTAR,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.concat_frames_in_dir_task import ConcatFramesInDirTask

STACK_UKBBPPP_RABGAP1L = ConcatFramesInDirTask.create(
    asset_id="ukbb_ppp_rabgap1l_sumstats_stacked",
    source_dir_task=UKBB_PPP_RABGAP1L_UNTAR,
    path_glob="*.gz",
    read_spec_for_frames=DataFrameReadSpec(DataFrameTextFormat(separator=" ")),
)
