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
    PppProteinFile,
)
from mecfs_bio.build_system.task.ppp_database.protein_sample_size_task import (
    DEFAULT_HEAD_BYTES,
    PppProteinRef,
    PppProteinSampleSizeTask,
)

_PPP_SUMSTATS_DIR = Path(ppp_assets.__file__).parent
# The two UKB-PPP cohorts (see regenerate_ppp_manifest.py). "Combined" is the full
# European cohort (discovery + replication, n = 52,363) and is the default; "discovery"
# is the smaller randomly selected baseline sub-cohort (n = 34,557), kept for the assets
# already built against those files.
EUR_COMBINED_PPP_MANIFEST_PATH = _PPP_SUMSTATS_DIR / "eur_combined_ppp_manifest.csv"
EUR_DISCOVERY_PPP_MANIFEST_PATH = _PPP_SUMSTATS_DIR / "eur_discovery_ppp_manifest.csv"

# Manifest columns (see regenerate_ppp_manifest.py). OID is the primary key: a few
# control proteins share a gene symbol across Olink panels, distinguished only by OID.
_GENE_COL = "gene"
_UNIPROT_COL = "uniprot"
_OID_COL = "oid"
_VERSION_COL = "version"
_PANEL_COL = "panel"
_SYNAPSE_ID_COL = "synapse_id"
_FILENAME_COL = "filename"


def _protein_file_from_row(row: dict) -> PppProteinFile:
    """Build the structured protein identity from a manifest row, checking that the derived
    tar filename matches the manifest's filename column (a guard against a malformed row)."""
    protein_file = PppProteinFile(
        gene=row[_GENE_COL],
        uniprot=row[_UNIPROT_COL],
        oid=row[_OID_COL],
        version=row[_VERSION_COL],
        panel=row[_PANEL_COL],
        synid=row[_SYNAPSE_ID_COL],
    )
    assert protein_file.tar_filename == row[_FILENAME_COL], (
        f"derived tar filename {protein_file.tar_filename} != manifest filename "
        f"{row[_FILENAME_COL]}"
    )
    return protein_file


@frozen
class PppSlimProteinTaskCollection:
    protein_tasks: tuple[BuildSlimProteinParquetTask, ...]

    def terminal_tasks(self) -> list[Task]:
        return list(self.protein_tasks)


def generate_ppp_slim_protein_tasks(
    index_task: Task,
    index_name: str,
    manifest_path: Path = EUR_COMBINED_PPP_MANIFEST_PATH,
    include_sample_size: bool = True,
) -> PppSlimProteinTaskCollection:
    """
    Build one BuildSlimProteinParquetTask per UKB-PPP protein listed in the manifest.

    index_task: the shared ConstructPppVariantIndexTask every protein aligns onto.
    index_name: a short label for the (index, cohort) combination (e.g. hapmap_3 or
        hapmap_3_eur_combined). It namespaces both the asset ids and the on-disk paths, so
        the same protein aligned onto a different index -- or drawn from a different cohort
        -- produces distinct, non-colliding assets.
    manifest_path: the committed per-cohort manifest; defaults to the full European
        (combined) cohort.
    include_sample_size: store the per-variant sample size N alongside beta/se (default
        on). Set False for cohorts whose slim files were already built without it, so those
        assets stay bit-identical to what is on disk.
    """
    manifest = pl.read_csv(manifest_path)
    protein_tasks = tuple(
        BuildSlimProteinParquetTask.create(
            index_task=index_task,
            protein=_protein_file_from_row(row),
            asset_id=f"ppp_slim_{index_name}_{row[_GENE_COL]}_{row[_OID_COL]}",
            index_name=index_name,
            include_sample_size=include_sample_size,
        )
        for row in manifest.iter_rows(named=True)
    )
    return PppSlimProteinTaskCollection(protein_tasks=protein_tasks)


def generate_ppp_sample_size_task(
    asset_id: str,
    manifest_path: Path = EUR_COMBINED_PPP_MANIFEST_PATH,
    head_bytes: int = DEFAULT_HEAD_BYTES,
) -> PppProteinSampleSizeTask:
    """Build the per-protein sample-size task from the manifest. The sample sizes are a
    property of the protein GWAS, independent of which variant index the database uses, so
    this takes no index. Defaults to the full European (combined) cohort manifest."""
    manifest = pl.read_csv(manifest_path)
    protein_refs = tuple(
        PppProteinRef(
            oid=row[_OID_COL], gene=row[_GENE_COL], synid=row[_SYNAPSE_ID_COL]
        )
        for row in manifest.iter_rows(named=True)
    )
    return PppProteinSampleSizeTask.create(
        asset_id=asset_id, protein_refs=protein_refs, head_bytes=head_bytes
    )
