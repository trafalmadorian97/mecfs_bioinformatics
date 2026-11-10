import shutil
import socket
from pathlib import Path

import polars as pl
import pysam
import pytest
from attrs import frozen

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.protein_coloc.csf_window_fetcher import (
    CSF_HARMONISED_COLUMNS,
    CsfWindowFetcher,
    HarmonisedLocation,
    UrlIndexFetcher,
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
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
)
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA

_ACCESSION = GcstAccession("GCST90421540")
_SEQ_ID = "2418-55"
_APTAMER = CsfAptamerFile(
    analyte=Analyte("X2418.55"),
    seq_id=SeqId(_SEQ_ID),
    accession=_ACCESSION,
    uniprot=UniProtId("P02649"),
    entrez_gene_symbol="APOE",
    entrez_gene_id=EntrezGeneId("348"),
)


@frozen(slots=True)
class _PathResolver:
    path: str

    def __call__(self, accession: GcstAccession) -> HarmonisedLocation:
        return HarmonisedLocation(data=self.path, index=f"{self.path}.tbi")


@frozen(slots=True)
class _CopyThenRemoveIndex:
    """Copies the index into scratch and removes the original, so a fetch only works if it
    reads the copy it was handed."""

    def __call__(self, index: str, dest_dir: Path) -> Path:
        dest = dest_dir / Path(index).name
        shutil.copy(index, dest)
        Path(index).unlink()
        return dest


def _row(pos: int, beta: str = "0.1", eaf: str = "0.3", n: str = "3400") -> list[str]:
    values = {
        "chromosome": "19",
        "base_pair_location": str(pos),
        "effect_allele": "A",
        "other_allele": "G",
        "beta": beta,
        "standard_error": "0.02",
        "effect_allele_frequency": eaf,
        "p_value": "1e-5",
        "variant_id": f"19_{pos}_G_A",
        "rsid": f"rs{pos}",
        "neg_log_10_p_value": "5",
        "rs_id": f"rs{pos}",
        "n": n,
        "hm_coordinate_conversion": "lo",
        "hm_code": "10",
    }
    return [values[c] for c in CSF_HARMONISED_COLUMNS]


def _tabix_fixture(tmp_path: Path, rows: list[list[str]]) -> str:
    plain = tmp_path / "study.h.tsv"
    lines = ["\t".join(CSF_HARMONISED_COLUMNS)] + ["\t".join(r) for r in rows]
    plain.write_text("\n".join(lines) + "\n")
    return pysam.tabix_index(
        str(plain), seq_col=0, start_col=1, end_col=1, line_skip=1, force=True
    )


def _fetch(source: str, window: LocusWindow, tmp_path: Path) -> pl.DataFrame:
    fetcher = CsfWindowFetcher(
        aptamers=(_APTAMER,),
        resolver=_PathResolver(source),
        index_fetcher=_CopyThenRemoveIndex(),
    )
    frames = fetcher.fetch(_SEQ_ID, [window], tmp_path, make_wf())
    return frames[window.locus_id]


def test_fetch_returns_window_rows_in_window_schema(tmp_path: Path) -> None:
    source = _tabix_fixture(
        tmp_path, [_row(100, n="3100"), _row(200, n="3200"), _row(300, n="3300")]
    )
    frame = _fetch(source, LocusWindow("L1", chrom=19, start=150, end=300), tmp_path)
    assert frame.schema == WINDOW_SCHEMA
    assert frame[GWASLAB_POS_COL].to_list() == [200, 300]
    assert frame[GWASLAB_SAMPLE_SIZE_COLUMN].to_list() == [3200.0, 3300.0]


def test_fetch_window_on_absent_contig_is_empty(tmp_path: Path) -> None:
    source = _tabix_fixture(tmp_path, [_row(100)])
    frame = _fetch(source, LocusWindow("LX", chrom=23, start=1, end=1000), tmp_path)
    assert frame.height == 0
    assert frame.schema == WINDOW_SCHEMA


def test_fetch_window_starting_before_chromosome_start(tmp_path: Path) -> None:
    # A trait locus near a chromosome end: lead minus the half-width is negative.
    source = _tabix_fixture(tmp_path, [_row(100), _row(200)])
    frame = _fetch(source, LocusWindow("L0", chrom=19, start=-400, end=150), tmp_path)
    assert frame[GWASLAB_POS_COL].to_list() == [100]


def test_fetch_na_fields_become_null(tmp_path: Path) -> None:
    source = _tabix_fixture(tmp_path, [_row(100, beta="NA", eaf="NA")])
    frame = _fetch(source, LocusWindow("L1", chrom=19, start=1, end=1000), tmp_path)
    assert frame[GWASLAB_BETA_COL].to_list() == [None]
    assert frame[GWASLAB_EFFECT_ALLELE_FREQ_COL].to_list() == [None]


def test_fetch_reads_the_index_it_downloaded_to_scratch(tmp_path: Path) -> None:
    source = _tabix_fixture(tmp_path, [_row(100)])
    frame = _fetch(source, LocusWindow("L1", chrom=19, start=1, end=1000), tmp_path)
    assert not Path(f"{source}.tbi").exists()
    assert frame.height == 1


def test_fetch_unreachable_source_raises_oserror(tmp_path: Path) -> None:
    missing = str(tmp_path / "missing.h.tsv.gz")
    with pytest.raises(OSError):
        _fetch(missing, LocusWindow("L1", chrom=19, start=1, end=1000), tmp_path)


def test_fetch_unknown_protein_id_asserts(tmp_path: Path) -> None:
    source = _tabix_fixture(tmp_path, [_row(100)])
    fetcher = CsfWindowFetcher(aptamers=(_APTAMER,), resolver=_PathResolver(source))
    with pytest.raises(AssertionError):
        fetcher.fetch(
            "9999-99",
            [LocusWindow("L1", chrom=19, start=1, end=10)],
            tmp_path,
            make_wf(),
        )


def test_fetch_region_pysam_rejects_is_not_retried_as_io_error(tmp_path: Path) -> None:
    # An end beyond what htslib can address is a programming error, not a read failure.
    source = _tabix_fixture(tmp_path, [_row(100)])
    window = LocusWindow("LBIG", chrom=19, start=1, end=2**63)
    with pytest.raises(ValueError):
        _fetch(source, window, tmp_path)


def test_index_download_from_unresponsive_server_times_out(tmp_path: Path) -> None:
    # The listening socket accepts the connection (kernel backlog) but never answers.
    with socket.socket() as server:
        server.bind(("127.0.0.1", 0))
        server.listen(1)
        port = server.getsockname()[1]
        with pytest.raises(OSError):
            UrlIndexFetcher(timeout_s=0.5)(
                f"http://127.0.0.1:{port}/study.h.tsv.gz.tbi", tmp_path
            )
