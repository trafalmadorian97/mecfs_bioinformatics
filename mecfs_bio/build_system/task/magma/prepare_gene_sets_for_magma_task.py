"""
Task to convert MSigDB gene sets to the MAGMA column-format annotation file.
"""

from pathlib import Path
from typing import Sequence

import pandas as pd
import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.magma_constants import MAGMA_GENE_COL
from mecfs_bio.constants.msigdb_columns import (
    EXACT_SOURCE,
    NCBI_IDS,
    STANDARD_NAME,
    SYSTEMATIC_NAME,
)
from mecfs_bio.constants.vocabulary_classes.gene_set import MSigDBGeneSetSpec
from mecfs_bio.util.gene_set.msigdb_lookup import _apply_spec_mask

logger = structlog.get_logger()


@frozen
class PrepareGeneSetsForMagmaTask(Task):
    """
    Reads MSigDB gene sets from a parquet database and writes a MAGMA
    column-format annotation file.

    The output is a tab-separated text file with a leading GENE column
    (Entrez IDs) followed by one binary (0/1) column per gene set.  This
    format is consumed by MagmaGeneSetAnalysisTask with set_or_covar="covar" (--gene-covar).

    All genes present anywhere in the MSigDB parquet are included as rows
    (required for competitive gene set analysis).  NCBI IDs that are None
    (unmapped genes) are silently dropped.
    """

    meta: Meta
    gene_sets: Sequence[MSigDBGeneSetSpec]
    parquet_db_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.parquet_db_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source = fetch(self.parquet_db_task.asset_id)
        assert isinstance(source, FileAsset)

        df = pd.read_parquet(
            source.path,
            columns=[STANDARD_NAME, SYSTEMATIC_NAME, EXACT_SOURCE, NCBI_IDS],
        )

        gs_to_ids: dict[str, frozenset[int]] = {}
        for spec in self.gene_sets:
            rows = df[_apply_spec_mask(df, spec)]
            if len(rows) != 1:
                raise ValueError(
                    f"Expected exactly 1 match for {STANDARD_NAME}={spec.standard_name!r}, "
                    f"got {len(rows)}"
                )
            ids = frozenset(
                int(x) for x in rows.iloc[0][NCBI_IDS] if x is not None and x == x
            )
            gs_to_ids[spec.standard_name] = ids

        all_genes = sorted(
            pl.from_pandas(df[[NCBI_IDS]])
            .select(pl.col(NCBI_IDS).explode().drop_nulls().cast(pl.Int64).unique())[
                NCBI_IDS
            ]
            .to_list()
        )
        logger.info(
            "building MAGMA gene set file",
            n_gene_sets=len(self.gene_sets),
            n_genes=len(all_genes),
        )

        col_data: dict[str, list] = {MAGMA_GENE_COL: all_genes}
        for spec in self.gene_sets:
            ids = gs_to_ids[spec.standard_name]
            col_data[spec.standard_name] = [1 if g in ids else 0 for g in all_genes]

        out_path = scratch_dir / "gene_sets_for_magma.txt"
        pd.DataFrame(col_data).to_csv(out_path, sep="\t", index=False)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        gene_sets: Sequence[MSigDBGeneSetSpec],
        parquet_db_task: Task,
    ) -> "PrepareGeneSetsForMagmaTask":
        source_meta = parquet_db_task.meta
        meta: Meta
        if isinstance(source_meta, ReferenceFileMeta):
            meta = ReferenceFileMeta(
                group=source_meta.group,
                sub_group=source_meta.sub_group,
                sub_folder=source_meta.sub_folder,
                id=AssetId(asset_id),
                extension=".txt",
            )
        elif isinstance(source_meta, SimpleFileMeta):
            meta = SimpleFileMeta(AssetId(asset_id))
        else:
            raise ValueError(
                f"Unsupported parquet_db_task meta type: {type(source_meta)}"
            )
        return cls(
            meta=meta,
            gene_sets=list(gene_sets),
            parquet_db_task=parquet_db_task,
        )
