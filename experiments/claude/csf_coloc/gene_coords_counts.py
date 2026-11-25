"""
One-off check of CsfGeneCoordsReader on the real manifest and the GRCh38 MAGMA gene table:
how many of the 7,008 aptamers get a gene interval (expected about 6,891 by id - 3 conflicts
+ 7 symbol fallbacks + the "1442|1442"-style repeated single ids).

Usage: pixi r python experiments/claude/csf_coloc/gene_coords_counts.py
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.asset_generator.csf_slim_aptamer_asset_generator import (
    load_csf_aptamer_files,
)
from mecfs_bio.assets.reference_data.magma_gene_locations.processed.magma_entrez_gene_locations_data_build_38_unzipped import (
    MAGMA_ENTREZ_GENE_LOCATION_REFERENCE_DATA_BUILD_38_EXTRACTED as MAGMA,
)
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.task.protein_coloc.csf_gene_coords_reader import (
    CsfGeneCoordsReader,
)


def main() -> None:
    aptamers = load_csf_aptamer_files()
    path = DEFAULT_RUNNER.meta_to_path(MAGMA.meta)
    table = (
        scan_dataframe_asset(FileAsset(path), meta=MAGMA.meta, parquet_backend="polars")
        .collect()
        .to_polars()
    )
    coords = CsfGeneCoordsReader(aptamers=aptamers)(table)
    print(f"aptamers {len(aptamers)}  with interval {len(coords)}")
    missing = [a for a in aptamers if a.seq_id not in coords]
    print(
        pl.DataFrame(
            {
                "seq_id": [a.seq_id for a in missing],
                "symbol": [a.entrez_gene_symbol for a in missing],
                "entrez_id": [a.entrez_gene_id for a in missing],
            }
        ).filter(pl.col("entrez_id").str.contains(r"^\d+$"))
    )


if __name__ == "__main__":
    main()
