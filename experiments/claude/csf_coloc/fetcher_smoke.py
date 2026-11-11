"""
Smoke test of CsfWindowFetcher against the real GWAS Catalog: fetch one APOE aptamer's
chr19:44.4-45.4 Mb window by remote tabix and print a summary.

Usage: pixi r python experiments/claude/csf_coloc/fetcher_smoke.py
"""

import tempfile
import time
from pathlib import Path

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.protein_coloc.csf_window_fetcher import (
    CsfWindowFetcher,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import LocusWindow
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.csf_database_constants import (
    Analyte,
    EntrezGeneId,
    GcstAccession,
    SeqId,
    UniProtId,
)

APOE = CsfAptamerFile(
    analyte=Analyte("X2418.55"),
    seq_id=SeqId("2418-55"),
    accession=GcstAccession("GCST90424891"),
    uniprot=UniProtId("P02649"),
    entrez_gene_symbol="APOE",
    entrez_gene_id=EntrezGeneId("348"),
)


def main() -> None:
    fetcher = CsfWindowFetcher(aptamers=(APOE,))
    window = LocusWindow("apoe", chrom=19, start=44_400_000, end=45_400_000)
    start = time.time()
    with tempfile.TemporaryDirectory() as tmp:
        frame = fetcher.fetch("2418-55", [window], Path(tmp), make_wf())["apoe"]
    print(f"fetched {frame.height} rows in {time.time() - start:.1f}s")
    print(frame.head())
    print(frame.null_count())
    print(frame.sort("BETA", descending=True).head(3))


if __name__ == "__main__":
    main()
