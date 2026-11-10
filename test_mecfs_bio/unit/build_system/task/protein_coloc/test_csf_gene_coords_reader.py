import polars as pl

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.protein_coloc.csf_gene_coords_reader import (
    CsfGeneCoordsReader,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_ENTREZ_ID_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.constants.csf_database_constants import (
    Analyte,
    EntrezGeneId,
    GcstAccession,
    SeqId,
    UniProtId,
)
from mecfs_bio.constants.protein_coloc_constants import ProteinId
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

# MAGMA-shaped rows: (entrez id, chrom, start, end, symbol), as the .loc file is read.
_MAGMA_ROWS = [
    (348, "19", 44_905_754, 44_909_393, "APOE"),
    (64115, "10", 71_747_556, 71_773_520, "VSIR"),
    (7322, "20", 50_081_124, 50_190_583, "TMEM189-UBE2V1"),
    (7335, "20", 50_100_000, 50_190_583, "UBE2V1"),
    (102724560, "21", 6_444_869, 6_467_509, "LOC102724560"),
    (875, "21", 43_053_191, 43_076_943, "CBS"),
    (2574, "X", 49_571_000, 49_582_000, "GAGE2B"),
]

_APOE = GenomicInterval(chrom=19, start=44_905_754, end=44_909_393)


def _magma() -> pl.DataFrame:
    return pl.DataFrame(
        {
            GENE_INFO_ENTREZ_ID_COL: [r[0] for r in _MAGMA_ROWS],
            GENE_INFO_CHROM_COL: [r[1] for r in _MAGMA_ROWS],
            GENE_INFO_START_COL: [r[2] for r in _MAGMA_ROWS],
            GENE_INFO_END_COL: [r[3] for r in _MAGMA_ROWS],
            GENE_INFO_STRAND_COL: ["+"] * len(_MAGMA_ROWS),
            GENE_INFO_NAME_COL: [r[4] for r in _MAGMA_ROWS],
        }
    )


def _aptamer(seq_id: str, entrez_id: str | None, symbol: str) -> CsfAptamerFile:
    return CsfAptamerFile(
        analyte=Analyte(f"X{seq_id.replace('-', '.')}"),
        seq_id=SeqId(seq_id),
        accession=GcstAccession("GCST90421540"),
        uniprot=UniProtId("P00000"),
        entrez_gene_symbol=symbol,
        entrez_gene_id=None if entrez_id is None else EntrezGeneId(entrez_id),
    )


def _coords(*aptamers: CsfAptamerFile) -> dict[ProteinId, GenomicInterval]:
    return CsfGeneCoordsReader(aptamers=aptamers)(_magma())


def test_entrez_id_found_gives_interval() -> None:
    # The aptamer table still calls VSIR by its old symbol C10orf54.
    coords = _coords(_aptamer("100-1", "64115", "C10orf54"))
    assert coords == {
        ProteinId("100-1"): GenomicInterval(chrom=10, start=71_747_556, end=71_773_520)
    }


def test_two_aptamers_same_gene_both_mapped() -> None:
    coords = _coords(
        _aptamer("2418-55", "348", "APOE"), _aptamer("2937-10", "348", "APOE")
    )
    assert coords == {ProteinId("2418-55"): _APOE, ProteinId("2937-10"): _APOE}


def test_id_and_symbol_nearby_keeps_id_interval() -> None:
    # Id points at the readthrough, symbol at the gene inside it: under 1 Mb apart.
    coords = _coords(_aptamer("200-1", "7322", "UBE2V1"))
    assert coords == {
        ProteinId("200-1"): GenomicInterval(chrom=20, start=50_081_124, end=50_190_583)
    }


def test_id_and_symbol_far_apart_has_no_interval() -> None:
    # Id points at the chr21 6.4 Mb duplicate copy, symbol at canonical CBS 36 Mb away.
    assert _coords(_aptamer("300-1", "102724560", "CBS")) == {}


def test_id_missing_falls_back_to_symbol() -> None:
    coords = _coords(_aptamer("400-1", "100101629", "GAGE2B"))
    assert coords == {
        ProteinId("400-1"): GenomicInterval(chrom=23, start=49_571_000, end=49_582_000)
    }


def test_multi_gene_id_has_no_interval() -> None:
    assert _coords(_aptamer("500-1", "348|875", "APOE|CBS")) == {}


def test_repeated_single_id_is_one_gene() -> None:
    assert _coords(_aptamer("510-1", "348|348", "APOE|APOE")) == {
        ProteinId("510-1"): _APOE
    }


def test_id_and_symbol_missing_from_magma_has_no_interval() -> None:
    assert _coords(_aptamer("600-1", "1", "NOTAGENE")) == {}
    assert _coords(_aptamer("700-1", None, "P12345")) == {}
