"""
DenseWindowFetcher for the Western et al. 2024 CSF pQTL GWAS: reads trait-locus windows of one
aptamer's GWAS Catalog harmonised summary statistics by tabix range requests, without
downloading the file.

Every study has a bgzipped, tabix-indexed GRCh38 harmonised file with a header line and the
columns in CSF_HARMONISED_COLUMNS. Contigs are "1".."22": there is no chrX, so a window on a
chromosome the index does not list yields an empty frame. Missing values are written "NA" and
become nulls. Effect allele and frequency are reported as the file has them; orientation to
the reference genome happens downstream, as for every source.

The tabix index is downloaded into the call's scratch directory and handed to pysam
explicitly; left to itself, htslib would save every remote index into the process's working
directory. Failures to download the index (including a stalled connection, after timeout_s)
or to open or read the source surface as OSError, so the coloc task's retry wrapper retries
them; a region pysam rejects is a programming error and stays a ValueError, failing at once.
"""

import shutil
import urllib.request
from collections.abc import Sequence
from pathlib import Path, PurePosixPath
from typing import Final, Protocol
from urllib.parse import urlparse

import polars as pl
import pysam
from attrs import frozen

from mecfs_bio.build_system.task.csf_database.build_slim_aptamer_parquet_task import (
    CsfAptamerFile,
)
from mecfs_bio.build_system.task.csf_database.gwas_catalog_url import (
    gwas_catalog_harmonised_url,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import LocusWindow
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.csf_database_constants import GcstAccession
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.protein_coloc_constants import WINDOW_SCHEMA

# Column contract of the Western et al. harmonised files (verified 2026-11-05).
CSF_HARMONISED_COLUMNS: Final = (
    "chromosome",
    "base_pair_location",
    "effect_allele",
    "other_allele",
    "beta",
    "standard_error",
    "effect_allele_frequency",
    "p_value",
    "variant_id",
    "rsid",
    "neg_log_10_p_value",
    "rs_id",
    "n",
    "hm_coordinate_conversion",
    "hm_code",
)
_HARMONISED_TO_WINDOW: Final = {
    "chromosome": GWASLAB_CHROM_COL,
    "base_pair_location": GWASLAB_POS_COL,
    "effect_allele": GWASLAB_EFFECT_ALLELE_COL,
    "other_allele": GWASLAB_NON_EFFECT_ALLELE_COL,
    "beta": GWASLAB_BETA_COL,
    "standard_error": GWASLAB_SE_COL,
    "effect_allele_frequency": GWASLAB_EFFECT_ALLELE_FREQ_COL,
    "n": GWASLAB_SAMPLE_SIZE_COLUMN,
}
_MISSING: Final = "NA"


@frozen(slots=True)
class HarmonisedLocation:
    """URL or local path of a bgzipped harmonised file (data) and of its tabix index."""

    data: str
    index: str


class HarmonisedLocationResolver(Protocol):
    def __call__(self, accession: GcstAccession) -> HarmonisedLocation: ...


@frozen(slots=True)
class GwasCatalogHarmonisedResolver:
    def __call__(self, accession: GcstAccession) -> HarmonisedLocation:
        url = gwas_catalog_harmonised_url(accession)
        return HarmonisedLocation(data=url, index=f"{url}.tbi")


class IndexFetcher(Protocol):
    def __call__(self, index: str, dest_dir: Path) -> Path:
        """Make the index at a URL or local path available as a local file in dest_dir."""
        ...


@frozen(slots=True)
class UrlIndexFetcher:
    """
    Downloads the index (a few MB) with urllib; URLError and socket timeouts are OSErrors.
    timeout_s: give up on a connection that sends nothing for this long.
    """

    timeout_s: float = 60.0

    def __call__(self, index: str, dest_dir: Path) -> Path:
        dest = dest_dir / PurePosixPath(urlparse(index).path).name
        with (
            urllib.request.urlopen(index, timeout=self.timeout_s) as response,
            open(dest, "wb") as out,
        ):
            shutil.copyfileobj(response, out)
        return dest


def _rows_to_window_frame(rows: list[list[str]]) -> pl.DataFrame:
    for row in rows:
        assert len(row) == len(CSF_HARMONISED_COLUMNS), (
            f"expected {len(CSF_HARMONISED_COLUMNS)} fields, got {len(row)}: {row}"
        )
    columns = {
        window_col: [None if row[i] == _MISSING else row[i] for row in rows]
        for i, name in enumerate(CSF_HARMONISED_COLUMNS)
        if (window_col := _HARMONISED_TO_WINDOW.get(name)) is not None
    }
    return (
        pl.DataFrame(columns, schema={c: pl.String for c in columns})
        .cast(WINDOW_SCHEMA)
        .select(WINDOW_SCHEMA.names())
    )


def _tabix_start(window: LocusWindow) -> int:
    """0-based start of the closed 1-based window, clipped at the chromosome start (a trait
    locus near a chromosome end has a negative start, which pysam rejects)."""
    return max(window.start - 1, 0)


def _read_window(tabix: pysam.TabixFile, window: LocusWindow) -> pl.DataFrame:
    if str(window.chrom) not in set(tabix.contigs):
        return _rows_to_window_frame([])
    # Building the iterator validates the region: a ValueError here is a bad request.
    lines = tabix.fetch(str(window.chrom), _tabix_start(window), window.end)
    try:
        return _rows_to_window_frame([line.split("\t") for line in lines])
    except ValueError as error:
        # pysam reports read failures (truncated stream, bad block) as ValueError.
        raise OSError(f"tabix read of {tabix.filename!r} failed: {error}") from error


def _open_tabix(source: str, index: Path) -> pysam.TabixFile:
    try:
        return pysam.TabixFile(source, index=str(index))
    except ValueError as error:
        raise OSError(f"opening {source} failed: {error}") from error


def _fetch_windows(
    source: str, index: Path, windows: Sequence[LocusWindow]
) -> dict[str, pl.DataFrame]:
    with _open_tabix(source, index) as tabix:
        return {window.locus_id: _read_window(tabix, window) for window in windows}


@frozen(slots=True)
class CsfWindowFetcher:
    """
    aptamers: the manifest rows of every aptamer that may be requested, looked up by seq id
        (the coloc protein id).
    resolver: where an aptamer's harmonised file lives; the GWAS Catalog FTP by default.
    index_fetcher: brings the tabix index into the scratch directory.
    """

    aptamers: tuple[CsfAptamerFile, ...]
    resolver: HarmonisedLocationResolver = GwasCatalogHarmonisedResolver()
    index_fetcher: IndexFetcher = UrlIndexFetcher()

    def _aptamer(self, protein_id: str) -> CsfAptamerFile:
        matches = [a for a in self.aptamers if a.seq_id == protein_id]
        assert len(matches) == 1, f"expected one manifest row for {protein_id}"
        return matches[0]

    def fetch(
        self,
        protein_id: str,
        windows: Sequence[LocusWindow],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[str, pl.DataFrame]:
        aptamer = self._aptamer(protein_id)
        location = self.resolver(aptamer.accession)
        index_dir = scratch_dir / f"tabix_index_{protein_id}"
        index_dir.mkdir(parents=True, exist_ok=True)
        index = self.index_fetcher(location.index, index_dir)
        try:
            return _fetch_windows(location.data, index, windows)
        finally:
            shutil.rmtree(index_dir)
