"""
Asset generator to automatically determine the optimal UKBB LD matrix to use for fine-mapping
"""

from pathlib import Path, PurePath

import pandas as pd
from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask

BROAD_UKBB_FILE_LIST = Path("mecfs_bio/vend_files/broad_ukbb_ld_matrix_file_list.txt")


def get_broad_ukbb_ld_matrix_file_info() -> pd.DataFrame:
    df = pd.read_csv(
        BROAD_UKBB_FILE_LIST, sep=r"\s+", names=[f"col_{i}" for i in range(6)]
    )
    split = df["col_3"].str.split(".", expand=True)
    split = split.iloc[2:-1]
    intervals = split[0].str.split("_", expand=True)
    intervals.columns = ["chrom", "start", "end"]
    intervals["chrom"] = intervals["chrom"].str.removeprefix("chr").astype(int)
    intervals["start"] = intervals["start"].astype(int)
    intervals["end"] = intervals["end"].astype(int)
    return intervals


@frozen
class GenomicInterval:
    chrom: int
    start: int
    end: int


def get_optimal_ukbb_ld_interval(chrom: int, pos: int) -> GenomicInterval:
    broad_ld_df = get_broad_ukbb_ld_matrix_file_info()

    valid_intervals_df = broad_ld_df.loc[
        (broad_ld_df["chrom"] == chrom)
        & (broad_ld_df["start"] <= pos)
        & (broad_ld_df["end"] >= pos)
    ]
    if len(valid_intervals_df) == 0:
        raise ValueError(
            f"Could not find a broad UKKBB LD matrix for chromosome {chrom}, position {pos}"
        )

    valid_intervals_df["start_dist"] = pos - valid_intervals_df["start"]
    valid_intervals_df["end_dist"] = valid_intervals_df["end"] - pos
    valid_intervals_df["centrality"] = (
        valid_intervals_df.loc[:, ["start_dist", "end_dist"]]
    ).min(axis=1)
    valid_intervals_df = valid_intervals_df.sort_values(
        by="centrality", ascending=False
    )
    optimal_interval = valid_intervals_df.iloc[0]
    return GenomicInterval(
        optimal_interval["chrom"].item(),
        optimal_interval["start"].item(),
        optimal_interval["end"].item(),
    )


def get_genomic_interval_stem_name(
    interval: GenomicInterval,
) -> str:
    return f"chr{interval.chrom}_{interval.start}_{interval.end}"


def get_genomic_interval_ld_labels_url(interval: GenomicInterval) -> str:
    return (
        "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/"
        + get_genomic_interval_stem_name(interval)
        + ".gz"
    )


def get_genomic_interval_ld_matrix_url(interval: GenomicInterval) -> str:
    return (
        "https://broad-alkesgroup-ukbb-ld.s3.amazonaws.com/UKBB_LD/"
        + get_genomic_interval_stem_name(interval)
        + ".npz"
    )


def get_genomic_interval_ld_labels_task(interval: GenomicInterval) -> Task:
    stem = get_genomic_interval_stem_name(interval)
    url = get_genomic_interval_ld_labels_url(interval)
    return DownloadFileTask(
        meta=ReferenceFileMeta(
            group="ukbb_reference_ld",
            sub_group=stem,
            sub_folder=PurePath("raw"),
            id=AssetId(f"ukbb_broad_ld_{stem}_labels"),
            filename=stem,
            extension=".gz",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator="\t")),
        ),
        url=url,
        md5_hash=None,
    )


def get_genomic_interval_ld_matrix_task(interval: GenomicInterval) -> Task:
    stem = get_genomic_interval_stem_name(interval)
    url = get_genomic_interval_ld_matrix_url(interval)
    return DownloadFileTask(
        meta=ReferenceFileMeta(
            group="ukbb_reference_ld",
            sub_group=stem,
            sub_folder=PurePath("raw"),
            id=AssetId(f"ukbb_broad_ld_{stem}_matrix"),
            filename=stem,
            extension=".npz",
        ),
        url=url,
        md5_hash=None,
    )


def get_ld_labels_and_matrix_task_for_chrom_pos_build_37(
    chrom: int, pos: int
) -> tuple[Task, Task]:
    interval = get_optimal_ukbb_ld_interval(chrom=chrom, pos=pos)
    return get_ld_labels_and_matrix_task_for_genomic_interval_build_37(
        interval=interval
    )


def get_ld_labels_and_matrix_task_for_genomic_interval_build_37(
    interval: GenomicInterval,
):
    return (
        get_genomic_interval_ld_labels_task(interval=interval),
        get_genomic_interval_ld_matrix_task(interval=interval),
    )
