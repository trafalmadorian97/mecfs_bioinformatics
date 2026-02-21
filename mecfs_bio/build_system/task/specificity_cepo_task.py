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
from mecfs_bio.build_system.task.specificity_frac_task import (
    filter_by_cell_count,
    filter_missing_genes,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class PrepareSpecificityCepo(Task):
    """
    Task to compute the specificity of genes for cell types using the CEPO specificity metric.

    The CEPO metric is based on differential stability: a gene is considered specific for a cell type if its expression in stable in that
    cell type, but not other cell types

    see
    Kim, Hani Jieun, et al. "Cepo uncovers cell identity through differential stability." bioRxiv (2021): 2021-01.
    """

    _meta: Meta
    long_count_df_task: Task
    cell_type_col: str
    count_col: str
    gene_col: str
    cell_col: str
    min_cells_per_type: int = 0
    epsilon: float = 0.0001
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
        num_genes = df.select(narwhals.col(self.gene_col).n_unique()).collect().item()
        df = _compute_inv_coef_prop_nonzero(
            df=df,
            cell_type_col=self.cell_type_col,
            count_col=self.count_col,
            gene_col=self.gene_col,
            epsilon=self.epsilon,
        )
        df = _compute_ranks(df=df, cell_type_col=self.cell_type_col)
        df = _compute_stability(
            df=df,
            num_genes=num_genes,
        )
        df = _compute_differential_stability(
            df=df,
            cell_type_col=self.cell_type_col,
            gene_col=self.gene_col,
        )
        df = df.select([self.cell_type_col, self.gene_col, DIFFERENTIAL_STABILITY])
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
        epsilon: float = 0.0001,
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
            epsilon=epsilon,
            out_format=out_format,
            pre_pipe=pre_pipe,
            post_pipe=post_pipe,
            cell_col=cell_col,
            min_cells_per_type=min_cells_per_type,
        )


_count_mean = "__count_mean__"
_count_std = "__count_std__"
_inv_coef = "__count_inv_coef_var__"
_count_prop_nonzero = "__count_prop_nonzero__"

_inv_coef_rank = "__inv_coef_rank__"
_nonzero_rank = "__nonzero_rank__"

_stability = "__stability__"
DIFFERENTIAL_STABILITY = "differential_stability"


def _compute_inv_coef_prop_nonzero(
    df: narwhals.LazyFrame,
    cell_type_col: str,
    count_col: str,
    gene_col: str,
    epsilon: float,
) -> narwhals.LazyFrame:
    df = df.group_by([cell_type_col, gene_col]).agg(
        narwhals.col(count_col).mean().alias(_count_mean),
        narwhals.col(count_col).std().alias(_count_std),
        (((narwhals.col(count_col) > 0).sum()) / (narwhals.col(count_col).len())).alias(
            _count_prop_nonzero
        ),
    )
    df = df.with_columns(
        (narwhals.col(_count_std) / (narwhals.col(_count_mean) + epsilon)).alias(
            _inv_coef
        )
    )
    return df


def _compute_ranks(
    df: narwhals.LazyFrame,
    cell_type_col: str,
) -> narwhals.LazyFrame:
    return df.with_columns(
        narwhals.col(_inv_coef)
        .rank("max")
        .over(
            cell_type_col,
            order_by=_inv_coef,
        )
        .alias(_inv_coef_rank),
        narwhals.col(_count_prop_nonzero)
        .rank("max")
        .over(cell_type_col, order_by=_count_prop_nonzero)
        .alias(_nonzero_rank),
    )


def _compute_stability(df: narwhals.LazyFrame, num_genes: int) -> narwhals.LazyFrame:
    return df.with_columns(
        (
            1
            - (narwhals.col(_inv_coef_rank) + narwhals.col(_count_prop_nonzero))
            / (2 * num_genes)
        ).alias(_stability)
    )


def _compute_differential_stability(
    df: narwhals.LazyFrame, cell_type_col: str, gene_col: str
) -> narwhals.LazyFrame:
    return df.with_columns(
        (
            narwhals.col(_stability) - narwhals.col(_stability).mean().over(gene_col)
        ).alias(DIFFERENTIAL_STABILITY)
    )
