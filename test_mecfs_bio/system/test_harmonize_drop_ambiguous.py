"""
System test of genome-reference harmonization on the first rows of DecodeME (lifted over
to build 37), using the real hg19 FASTA and 1000 Genomes EUR panel assets.
"""

import tempfile
from pathlib import Path

import polars as pl

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.filtered_snps_gwas_1 import (
    DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.thousand_genomes.eur_panel_allele_frequencies import (
    THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.dataframe_output import ParquetOutFormat
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.head_pipe import HeadPipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
)

_decode_me_first_rows = PipeDataFrameTask.create(
    source_task=DECODE_ME_FILTER_SNPS_GWAS_1_TASK,
    asset_id="testing_decode_me_slice_task",
    out_format=ParquetOutFormat(),
    pipes=[HeadPipe(num_rows=10_000)],
)

_decode_me_first_rows_liftover_to_37 = GWASLabCreateSumstatsTask(
    df_source_task=_decode_me_first_rows,
    target_asset_id=AssetId(
        "testing_first_rows_decode_me_gwas_1_sumstats_liftover_to_37"
    ),
    basic_check=True,
    genome_build="infer",
    liftover_to="19",
)

_pre_harmonization_table = GwasLabSumstatsToTableTask.create_from_source_task(
    source_tsk=_decode_me_first_rows_liftover_to_37,
    asset_id="testing_first_rows_pre_harmonization_table",
    sub_dir="processed",
)

_harmonized_task = GenomeReferenceHarmonizationTask.create(
    asset_id="testing_first_rows_genome_reference_harmonized",
    sumstats_task=_pre_harmonization_table,
    fasta_task=UCSC_HG19_INDEXED_FASTA,
    panel_task=THOUSAND_GENOMES_EUR_HG19_PANEL_ALLELE_FREQUENCIES,
)


def test_genome_reference_harmonization_orients_every_kept_variant_to_the_reference():
    with tempfile.TemporaryDirectory() as tempdirname:
        tempdir = Path(tempdirname)
        asset_root = tempdir / "asset_store"
        asset_root.mkdir(parents=True, exist_ok=True)
        runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=tempdir / "info_store.yaml",
            asset_root=asset_root,
        )
        assets = runner.run([_harmonized_task, UCSC_HG19_INDEXED_FASTA])
        harmonized = assets[_harmonized_task.asset_id]
        fasta_asset = assets[UCSC_HG19_INDEXED_FASTA.asset_id]
        assert isinstance(harmonized, FileAsset)
        assert isinstance(fasta_asset, DirectoryAsset)
        table = pl.read_parquet(harmonized.path)
        assert table.height > 0
        fasta = IndexedFasta.open(fasta_asset.path)
        for chrom, group in table.group_by(GWASLAB_CHROM_COL):
            matches = reference_matches(
                fasta,
                chrom=int(chrom[0]),
                positions=group[GWASLAB_POS_COL].to_numpy(),
                alleles=group[GWASLAB_NON_EFFECT_ALLELE_COL],
            )
            assert matches.all()
