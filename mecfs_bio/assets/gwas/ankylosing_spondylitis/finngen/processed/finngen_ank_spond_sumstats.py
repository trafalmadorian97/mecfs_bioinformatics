from mecfs_bio.assets.gwas.ankylosing_spondylitis.finngen.raw.raw_finngen_ank_spond_data import \
    FINNGEN_ANKYLOSING_SPONDYLITIS_DATA_RAW
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import GWASLabCreateSumstatsTask, \
    GWASLabColumnSpecifiers
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.str_split_exact_col import SplitExactColPipe

FINGNEN_ANK_SPOND_SUMSTATS = GWASLabCreateSumstatsTask(
    df_source_task=FINNGEN_ANKYLOSING_SPONDYLITIS_DATA_RAW,
    asset_id=AssetId("finngen_ank_spond_sumstats"),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
    fmt=GWASLabColumnSpecifiers(
rsid="ID",
        chrom="CHROM",
        pos="POS",
        nea="REF",
        ea="ALT",
        beta="ES",
        se="SE",
        mlog10p="LP",
        eaf="AF"
    ),
pre_pipe=CompositePipe(
    [
        SplitExactColPipe(
            col_to_split="finn-b-M13_ANKYLOSPON",
            split_by=":",
            new_col_names=tuple(["ES","SE","LP","AF"])
        )
    ]
)
    # pre_pipe=pre_pipe,
)
