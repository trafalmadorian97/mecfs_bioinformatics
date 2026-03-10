import numpy as np
import pandas as pd
from pathlib import Path

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import DataFrameReadSpec, DataFrameTextFormat
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import ReferenceFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class FuseLociForLavaTask(Task):
    _meta: Meta
    original_loci_task: Task
    original_loci_per_fused_locus:int

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.original_loci_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_loci_asset = fetch(self.original_loci_task.asset_id)
        df =scan_dataframe_asset(source_loci_asset, meta=self.original_loci_task.meta).collect().to_pandas()
        n= self.original_loci_per_fused_locus
        out_frames= []
        for chr,grp in df.groupby("CHR"):
            chrom_rows = len(grp)
            for i in range(0, chrom_rows + n, n):
                to_fuse = grp.iloc[i:i + n]
                if len(to_fuse) > 0:
                    out_frames.append(
                        pd.DataFrame(
                            {
                                "CHR": [chr],
                                "START": [to_fuse["START"].min()],
                                "STOP": [to_fuse["STOP"].max()]
                            }
                        )
                    )
        result: pd.DataFrame = pd.concat(out_frames, ignore_index=True)
        result["LOC"] = (np.arange(len(result))+1)
        result  = result[["LOC","CHR","START","STOP"]]
        out_path = scratch_dir/"result.locfile"
        result.to_csv(out_path, index=False, sep=" ")
        return FileAsset(out_path)

    @classmethod
    def create(cls,
asset_id: str,
               original_loci_task: Task,
        original_loci_per_fused_locus: int,
    ):
        source_meta=original_loci_task.meta
        assert isinstance(source_meta, ReferenceFileMeta)
        meta= ReferenceFileMeta(
            group=source_meta.group,
            sub_group=source_meta.sub_group,
            sub_folder=source_meta.sub_folder,
            extension=".locfile",
            id=AssetId(asset_id),
            read_spec=DataFrameReadSpec(
                format=DataFrameTextFormat(separator=" ")
            )
        )
        return cls(
            meta=meta,
            original_loci_task=original_loci_task,
            original_loci_per_fused_locus=original_loci_per_fused_locus
        )

