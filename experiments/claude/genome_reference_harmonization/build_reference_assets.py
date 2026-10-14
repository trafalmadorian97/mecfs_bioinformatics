"""Materialize the genome-reference harmonization reference assets and report their size.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.build_reference_assets \
    2>&1 | tee experiments/claude/genome_reference_harmonization/build_reference_assets.log
"""

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg38_fasta import (
    UCSC_HG38_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL

FASTAS = [UCSC_HG19_INDEXED_FASTA, UCSC_HG38_INDEXED_FASTA]
PANELS = [
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
    THOUSAND_GENOMES_EUR_HG38_PANEL_ALLELE_FREQUENCIES,
]


def main() -> None:
    assets = DEFAULT_RUNNER.run([*FASTAS, *PANELS])
    for task in FASTAS:
        asset = assets[task.asset_id]
        assert isinstance(asset, DirectoryAsset)
        fasta = IndexedFasta.open(asset.path)
        print(task.asset_id, "chromosomes:", sorted(fasta.entries))
    for task in PANELS:
        asset = assets[task.asset_id]
        assert isinstance(asset, FileAsset)
        per_chrom = (
            pl.scan_parquet(asset.path)
            .group_by(GWASLAB_CHROM_COL)
            .len()
            .sort(GWASLAB_CHROM_COL)
            .collect()
        )
        print(task.asset_id, f"{asset.path.stat().st_size / 1e9:.2f} GB")
        print(per_chrom)


if __name__ == "__main__":
    main()
