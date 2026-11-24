import polars as pl

from mecfs_bio.build_system.task.decode_database.decode_somascan_protein import (
    DecodeSomascanProtein,
)
from mecfs_bio.build_system.task.protein_coloc.decode_gene_coords_reader import (
    DecodeGeneCoordsReader,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_ENTREZ_ID_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.s3_access import S3Endpoint, S3ObjectRef
from mecfs_bio.constants.csf_database_constants import EntrezGeneId, SeqId, UniProtId
from mecfs_bio.constants.protein_coloc_constants import ProteinId
from mecfs_bio.constants.vocabulary_classes.genomic_interval import GenomicInterval

# (entrez id, chrom, start, end, symbol), as the MAGMA .loc file is read.
_MAGMA_ROWS = [
    (348, "19", 44_905_754, 44_909_393, "APOE"),
    (102724560, "21", 6_444_869, 6_467_509, "LOC102724560"),
    (875, "21", 43_053_191, 43_076_943, "CBS"),
    (2574, "X", 49_571_000, 49_582_000, "GAGE2B"),
]
_ENDPOINT = S3Endpoint(url="https://fake.example", region="r", credentials_profile="p")


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


def _protein(seq_id: str, entrez_id: str | None, symbol: str) -> DecodeSomascanProtein:
    return DecodeSomascanProtein(
        seq_id=SeqId(seq_id),
        gene_symbol=symbol,
        target_name=symbol,
        uniprot=UniProtId("P00000"),
        entrez_gene_id=None if entrez_id is None else EntrezGeneId(entrez_id),
        ref=S3ObjectRef(endpoint=_ENDPOINT, bucket="b", key=f"{seq_id}.txt.gz"),
        size=1,
        md5="0" * 32,
    )


def test_decode_gene_coords_uses_entrez_then_symbol_rules() -> None:
    proteins = (
        _protein("1-1", "348", "APOE"),
        # Id points at the chr21 duplicate copy, symbol at canonical CBS 36 Mb away.
        _protein("2-1", "102724560", "CBS"),
        _protein("3-1", "348|875", "APOE"),
        # Retired id: the symbol names one MAGMA gene.
        _protein("4-1", "100101629", "GAGE2B"),
    )
    assert DecodeGeneCoordsReader(proteins=proteins)(_magma()) == {
        ProteinId("1-1"): GenomicInterval(chrom=19, start=44_905_754, end=44_909_393),
        ProteinId("4-1"): GenomicInterval(chrom=23, start=49_571_000, end=49_582_000),
    }
