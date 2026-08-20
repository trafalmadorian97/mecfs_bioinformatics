"""
Asset generator for the compact Western et al. 2024 CSF pQTL store: one slim
per-aptamer aligned beta/se/N task for every aptamer in the manifest.

The CSF analogue of ukbb_ppp_slim_protein_asset_generator. Asset ids and paths are
namespaced by the SomaScan seq id (the aptamer primary key) rather than the gene
symbol, because gene symbol is not unique across aptamers.
"""

from pathlib import Path

import polars as pl
from attrs import frozen

import mecfs_bio.assets.reference_data.csf_pqtl_sumstats as csf_assets
from mecfs_bio.assets.reference_data.csf_pqtl_sumstats.regenerate_csf_manifest import (
    ACCESSION_MANIFEST_COLUMN,
    ANALYTE_MANIFEST_COLUMN,
    ENTREZ_SYMBOL_MANIFEST_COLUMN,
    MD5_MANIFEST_COLUMN,
    SEQ_ID_MANIFEST_COLUMN,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    BuildSlimCsfAptamerParquetTask,
    CsfAptamerFile,
)
from mecfs_bio.constants.csf_database_constants import (
    Analyte,
    GcstAccession,
    SeqId,
)

_CSF_SUMSTATS_DIR = Path(csf_assets.__file__).parent
CSF_APTAMER_MANIFEST_PATH = _CSF_SUMSTATS_DIR / "csf_aptamer_manifest.csv"


def _aptamer_file_from_row(row: dict) -> CsfAptamerFile:
    return CsfAptamerFile(
        analyte=Analyte(row[ANALYTE_MANIFEST_COLUMN]),
        seq_id=SeqId(row[SEQ_ID_MANIFEST_COLUMN]),
        accession=GcstAccession(row[ACCESSION_MANIFEST_COLUMN]),
        entrez_gene_symbol=row[ENTREZ_SYMBOL_MANIFEST_COLUMN],
    )


# The seq id (dash form, e.g. 13681-173) is not a safe asset-id token: '-' collides
# with the id delimiter. The analyte's dot form (X13681.173) is not either. Use the
# seq id with its separators normalized to underscores, which is unique per aptamer.
def _asset_id_token(seq_id: str) -> str:
    return seq_id.replace("-", "_")


@frozen
class CsfSlimAptamerTaskCollection:
    aptamer_tasks: tuple[BuildSlimCsfAptamerParquetTask, ...]

    def terminal_tasks(self) -> list[Task]:
        return list(self.aptamer_tasks)


def generate_csf_slim_aptamer_tasks(
    index_task: Task,
    index_name: str,
    manifest_path: Path = CSF_APTAMER_MANIFEST_PATH,
) -> CsfSlimAptamerTaskCollection:
    """
    Build one BuildSlimCsfAptamerParquetTask per aptamer listed in the manifest.

    index_task: the shared ConstructCsfVariantIndexTask every aptamer aligns onto.
    index_name: a short label for the index (e.g. hapmap_3). It namespaces both the
        asset ids and the on-disk paths, so the same aptamer aligned onto a different
        index produces distinct, non-colliding assets.
    manifest_path: the committed aptamer manifest.
    """
    manifest = pl.read_csv(manifest_path)
    aptamer_tasks = tuple(
        BuildSlimCsfAptamerParquetTask.create(
            index_task=index_task,
            aptamer=_aptamer_file_from_row(row),
            asset_id=(
                f"csf_slim_{index_name}_{row[ENTREZ_SYMBOL_MANIFEST_COLUMN]}_"
                f"{_asset_id_token(row[SEQ_ID_MANIFEST_COLUMN])}"
            ),
            index_name=index_name,
            md5_hash=row[MD5_MANIFEST_COLUMN],
        )
        for row in manifest.iter_rows(named=True)
    )
    return CsfSlimAptamerTaskCollection(aptamer_tasks=aptamer_tasks)
