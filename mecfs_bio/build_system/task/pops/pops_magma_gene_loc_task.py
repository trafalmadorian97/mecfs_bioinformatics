"""
Task to derive a MAGMA gene-location file from POPs' gene_annot_jun10.txt.

POPs requires that the MAGMA gene set be a subset of its own gene annotation
(pops.py indexes the gene annotation by every MAGMA gene). Running MAGMA with our
standard Ensembl v102 gene locations produces genes POPs does not know about, so a
POPs-consumed MAGMA analysis must instead use gene locations derived from POPs'
gene_annot_jun10.txt.

That annotation ships inside the POPs source tarball as a tab-separated table with
ENSGID, NAME, CHR, START, END, and TSS columns. This task rewrites it into MAGMA's
gene-location format (ENSGID, CHR, START, END, STRAND, NAME), inferring strand from
whether the TSS coincides with the start (+) or end (-) of the gene.
"""

from pathlib import Path

import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pops.pops_utils import GENE_ANNOT_RELATIVE_PATH
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()

# Columns of POPs' gene_annot_jun10.txt (tab-separated, with header).
_ANNOT_HEADER = ["ENSGID", "NAME", "CHR", "START", "END", "TSS"]


@frozen
class PopsMagmaGeneLocTask(Task):
    """Convert POPs' gene_annot_jun10.txt into a MAGMA gene-location file."""

    meta: Meta
    pops_source_task: Task

    @property
    def deps(self) -> list["Task"]:
        return [self.pops_source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.pops_source_task.asset_id)
        assert isinstance(source_asset, DirectoryAsset)
        annot_path = source_asset.path / GENE_ANNOT_RELATIVE_PATH
        out_path = scratch_dir / "pops_gene_annot.gene.loc"
        num_genes = write_magma_gene_loc(annot_path, out_path)
        logger.info(
            "Wrote MAGMA gene-location file from POPs annotation", genes=num_genes
        )
        return FileAsset(out_path)


def write_magma_gene_loc(annot_path: Path, out_path: Path) -> int:
    """Convert a POPs gene_annot_jun10.txt table into a MAGMA gene-location file.

    Writes tab-separated ENSGID, CHR, START, END, STRAND, NAME rows (no header),
    matching the format MAGMA's --gene-loc expects. Returns the number of genes
    written. Strand is + when the TSS is the start coordinate, - otherwise.
    """
    with (
        annot_path.open("r", encoding="utf-8") as source,
        out_path.open("w", encoding="utf-8") as out_file,
    ):
        header = source.readline().rstrip("\n").split("\t")
        assert header == _ANNOT_HEADER, f"Unexpected annotation header: {header}"
        num_genes = 0
        for line in source:
            ensgid, name, chrom, start, end, tss = line.rstrip("\n").split("\t")
            strand = "+" if tss == start else "-"
            out_file.write("\t".join([ensgid, chrom, start, end, strand, name]) + "\n")
            num_genes += 1
    return num_genes
