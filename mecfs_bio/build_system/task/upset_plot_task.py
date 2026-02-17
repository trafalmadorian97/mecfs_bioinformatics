from pathlib import Path, PurePath
from typing import Sequence

import matplotlib.pyplot as plt
from attrs import frozen
from upsetplot import UpSet, from_contents

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.plot_file_meta import GWASPlotFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec
from mecfs_bio.build_system.meta.read_spec.read_dataframe import (
    scan_dataframe,
    scan_dataframe_asset,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


@frozen
class FileSetSource:
    name: str
    task: Task
    col_name: str
    pipe: DataProcessingPipe = IdentityPipe()


@frozen
class DirSetSource:
    name: str
    task: Task
    file_in_dir: PurePath
    read_spec: DataFrameReadSpec
    col_name: str
    pipe: DataProcessingPipe = IdentityPipe()


SetSource = FileSetSource | DirSetSource


def load_contents(set_source: SetSource, fetch: Fetch) -> list[str]:
    if isinstance(set_source, FileSetSource):
        asset = fetch(set_source.task.asset_id)
        df = (
            set_source.pipe.process(
                scan_dataframe_asset(asset, meta=set_source.task.meta)
            )
            .collect()
            .to_pandas()
        )
        return df[set_source.col_name].tolist()
    elif isinstance(set_source, DirSetSource):
        asset = fetch(set_source.task.asset_id)
        assert isinstance(asset, DirectoryAsset)
        df = (
            set_source.pipe.process(
                scan_dataframe(
                    path=asset.path / set_source.file_in_dir, spec=set_source.read_spec
                )
            )
            .collect()
            .to_pandas()
        )
        return df[set_source.col_name].tolist()
    else:
        raise ValueError("Unknown set source")


@frozen
class UpSetPlotTask(Task):
    _meta: Meta
    set_sources: Sequence[SetSource]

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result = [item.task for item in self.set_sources]
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        contents_dict = {
            item.name: load_contents(item, fetch=fetch) for item in self.set_sources
        }
        sets = from_contents(contents_dict)
        # fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(16, 10))
        UpSet(
            sets,
            show_counts=True,
        ).plot(
            # fig
        )
        write_plots_to_dir(
            scratch_dir,
            {
                "upset": plt.gcf(),
            },
        )
        return FileAsset(scratch_dir / "upset.png")

    @classmethod
    def create(cls, asset_id: str, set_sources: Sequence[SetSource]):
        assert len(set_sources) >= 1
        source_meta = set_sources[0].task.meta
        if isinstance(source_meta, ResultDirectoryMeta):
            meta = GWASPlotFileMeta(
                trait=source_meta.trait,
                project=source_meta.project,
                extension=".png",
                id=AssetId(asset_id),
            )
            return cls(
                meta=meta,
                set_sources=set_sources,
            )
        raise ValueError(f"Unknown source meta {source_meta}")
