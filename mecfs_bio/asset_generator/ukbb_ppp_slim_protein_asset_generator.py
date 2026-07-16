"""
Asset generator for the compact UKB-PPP store: one slim per-protein aligned
beta/se task for every protein in the manifest.

"""

from pathlib import Path

import polars as pl
from attrs import frozen

import mecfs_bio.assets.reference_data.ukbb_ppp_sumstats as ppp_assets
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.ppp_database.build_slim_protein_parquet_task import (
    BuildSlimProteinParquetTask,
)

PPP_MANIFEST_PATH = Path(ppp_assets.__file__).parent / "ppp_manifest.csv"

# Manifest columns (see regenerate_ppp_manifest.py). OID is the primary key: a few
# control proteins share a gene symbol across Olink panels, distinguished only by OID.
_GENE_COL = "gene"
_OID_COL = "oid"
_SYNAPSE_ID_COL = "synapse_id"
_FILENAME_COL = "filename"


@frozen
class PppSlimProteinTaskCollection:
    protein_tasks: tuple[BuildSlimProteinParquetTask, ...]

    def terminal_tasks(self) -> list[Task]:
        return list(self.protein_tasks)


def generate_ppp_slim_protein_tasks(
    index_task: Task,
    index_name: str,
    manifest_path: Path = PPP_MANIFEST_PATH,
) -> PppSlimProteinTaskCollection:
    """
    Build one BuildSlimProteinParquetTask per UKB-PPP protein listed in the manifest.

    index_task: the shared ConstructPppVariantIndexTask every protein aligns onto.
    index_name: a short label for the index (e.g. hapmap3 or common_1kg). It namespaces
        both the asset ids and the on-disk paths, so the same protein aligned onto two
        different indices produces distinct, non-colliding assets.
    manifest_path: the committed ppp_manifest.csv; defaults to the packaged manifest.
    """
    manifest = pl.read_csv(manifest_path)
    protein_tasks = tuple(
        BuildSlimProteinParquetTask.create(
            index_task=index_task,
            synid=row[_SYNAPSE_ID_COL],
            expected_tar_filename=row[_FILENAME_COL],
            gene=row[_GENE_COL],
            asset_id=f"ppp_slim_{index_name}_{row[_GENE_COL]}_{row[_OID_COL]}",
            index_name=index_name,
        )
        for row in manifest.iter_rows(named=True)
    )
    return PppSlimProteinTaskCollection(protein_tasks=protein_tasks)
