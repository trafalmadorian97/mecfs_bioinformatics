"""Reference asset: the Eagle hg19 (build 37) genetic map and its parsed parquet.

Provides the per-position recombination rate (cM/Mb) used by the polyfun
explainability plot's recombination track. hg19 matches the polyfun LD panel.
"""

from pathlib import PurePath

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.download_file_task import DownloadFileTask
from mecfs_bio.build_system.task.genetic_map.parse_genetic_map_task import (
    ParseHg19GeneticMapTask,
)

GENETIC_MAP_HG19_RAW = DownloadFileTask(
    meta=ReferenceFileMeta(
        group="genetic_map",
        sub_group="hg19",
        sub_folder=PurePath("raw"),
        id=AssetId("genetic_map_hg19_eagle"),
        extension=".txt.gz",
    ),
    url="https://storage.googleapis.com/broad-alkesgroup-public/Eagle/downloads/tables/genetic_map_hg19_withX.txt.gz",
    md5_hash="930ba8e1435d54f68fb7a723fd3f0fa4",
)

GENETIC_MAP_HG19 = ParseHg19GeneticMapTask.create(
    asset_id="genetic_map_hg19",
    raw_map_task=GENETIC_MAP_HG19_RAW,
)
