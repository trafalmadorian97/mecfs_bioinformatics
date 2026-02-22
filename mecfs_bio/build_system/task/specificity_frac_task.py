from pathlib import Path

import narwhals
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipe_dataframe_task import (
    CSVOutFormat,
    OutFormat,
    ParquetOutFormat,
    get_extension_and_read_spec_from_format,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class PrepareSpecificityFraction(Task):
    """
    Task to compute the specificity of genes for cell types using the fractional specificity metric.

    In this metric, the specificity of a gene for a cell type is
    (mean expression in cell type)/(sum over all cell types of mean expression in those cell types)
    """

    _meta: Meta
    long_count_df_task: Task
    cell_type_col: str
    count_col: str
    gene_col: str
    cell_col: str
    min_cells_per_type: int = 0
    out_format: OutFormat = ParquetOutFormat()
    pre_pipe: DataProcessingPipe = IdentityPipe()
    post_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.long_count_df_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        long_count_df_asset = fetch(self.long_count_df_task.asset_id)

        df = scan_dataframe_asset(
            long_count_df_asset, meta=self.long_count_df_task.meta
        )
        df = self.pre_pipe.process(df)
        df = filter_by_cell_count(
            df,
            cell_type_col=self.cell_type_col,
            cell_col=self.cell_col,
            min_cells=self.min_cells_per_type,
        )
        df = filter_missing_genes(df, gene_col=self.gene_col, count_col=self.count_col)
        df = _compute_cell_type_means(
            df,
            cell_type_col=self.cell_type_col,
            count_col=self.count_col,
            gene_col=self.gene_col,
        )
        df = _compute_normalized_mean(
            df=df,
            gene_col=self.gene_col,
        )
        df = df.select([self.cell_type_col, self.gene_col, NORMALIZED_MEAN])
        df = self.post_pipe.process(df)
        out_path = scratch_dir / "out"
        if isinstance(self.out_format, CSVOutFormat):
            df.collect().to_pandas().to_csv(
                out_path, index=False, sep=self.out_format.sep
            )
        elif isinstance(self.out_format, ParquetOutFormat):
            df.sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        long_count_df_task: Task,
        cell_type_col: str,
        count_col: str,
        gene_col: str,
        cell_col: str,
        min_cells_per_type: int,
        out_format: OutFormat = ParquetOutFormat(),
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        post_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        extension, read_spec = get_extension_and_read_spec_from_format(
            out_format=out_format
        )
        source_meta = long_count_df_task.meta
        if isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=source_meta.sub_folder,
                id=AssetId(asset_id),
                filename=None,
                extension=extension,
                read_spec=read_spec,
            )
        else:
            raise ValueError(f"Unknown meta: {source_meta}")
        return cls(
            meta=meta,
            long_count_df_task=long_count_df_task,
            cell_type_col=cell_type_col,
            count_col=count_col,
            gene_col=gene_col,
            out_format=out_format,
            pre_pipe=pre_pipe,
            post_pipe=post_pipe,
            cell_col=cell_col,
            min_cells_per_type=min_cells_per_type,
        )


_mean_col = "__mean__"
NORMALIZED_MEAN = "normalized_mean"


def filter_missing_genes(
    df: narwhals.LazyFrame,
    gene_col: str,
    count_col: str,
) -> narwhals.LazyFrame:
    nz = df.group_by(gene_col).agg(
        (narwhals.col(count_col) > 0).sum().alias("nonzero_count")
    )
    nz = nz.filter(narwhals.col("nonzero_count") > 0)
    return df.join(nz, on=gene_col)


def filter_by_cell_count(
    df: narwhals.LazyFrame, cell_type_col: str, cell_col: str, min_cells: int
) -> narwhals.LazyFrame:
    cc = df.group_by(cell_type_col).agg(
        narwhals.col(cell_col).n_unique().alias("__cell_count")
    )
    cc = cc.filter(narwhals.col("__cell_count") >= min_cells)
    return df.join(cc, on=cell_type_col)


def _compute_cell_type_means(
    df: narwhals.LazyFrame,
    cell_type_col: str,
    gene_col: str,
    count_col: str,
) -> narwhals.LazyFrame:
    return df.group_by(
        [cell_type_col, gene_col],
    ).agg(narwhals.col(count_col).mean().alias(_mean_col))


def _compute_normalized_mean(
    df: narwhals.LazyFrame,
    gene_col: str,
) -> narwhals.LazyFrame:
    return df.with_columns(
        (
            narwhals.col(_mean_col) / (narwhals.col(_mean_col).sum().over(gene_col))
        ).alias(NORMALIZED_MEAN)
    )
