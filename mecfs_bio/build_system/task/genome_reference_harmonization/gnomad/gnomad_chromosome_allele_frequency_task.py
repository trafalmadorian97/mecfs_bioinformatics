"""
One chromosome of a gnomAD sites VCF, reduced to an allele-frequency table.

Streams the release's per-chromosome VCF (usually over HTTPS) with a single bcftools query
process, keeping PASS records polymorphic in at least one main group, and writes CHR, POS,
REF, ALT and AF_g, AN_g for every group of the release. One Task per chromosome keeps the
genome-wide build resumable: a failure costs one chromosome.

The header is checked first (every AF_g and AN_g declared, every contig line carrying the
release's assembly, this chromosome's contig declared), so a wrong release description fails
in seconds rather than after a long stream. bcftools query runs as one process rather than
a view | query pipe: execute_command runs commands through sh, where a pipeline's exit status
is that of its last command, so an upstream network failure could otherwise leave a
truncated table and exit 0. A single process lets htslib's BGZF integrity checks fail the
Task. Each retry overwrites the scratch TSV.

The TSV becomes parquet in bounded memory through the shared panel batch checks. For every
group, AF must be null exactly where AN is 0. Any REF that differs from a pure A/C/G/T FASTA
span is fatal; records over N or other IUPAC bases are dropped and counted. AN columns are
written with byte-stream-split encoding and AF columns with dictionary encoding, the smallest
combination measured.
"""

import re
import shlex
from collections.abc import Iterator
from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    gwaslab_code_to_contig_name,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    load_fasta,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.gnomad.gnomad_release import (
    GNOMAD_GROUP,
    GnomadRelease,
    gnomad_contig_name,
    gnomad_vcf_url,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.panel_batch_checks import (
    DEFAULT_BATCH_ROWS,
    SOURCE_CONTIG_COL,
    PanelBatchContext,
    PanelWriteSummary,
    write_checked_panel,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.allele_frequency_panel_constants import (
    GnomadGroup,
    panel_af_col,
    panel_an_col,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_POS_COL
from mecfs_bio.util.subproc.run_command import execute_command_with_retries

logger = structlog.get_logger()

DEFAULT_MAX_ATTEMPTS = 6
GNOMAD_PART_SUB_FOLDER = "per_chromosome"
_TSV_NULL = "."
_INFO_ID = re.compile(r"^##INFO=<ID=([^,>]+)")
_CONTIG_ID = re.compile(r"^##contig=<ID=([^,>]+)")
_ASSEMBLY = re.compile(r"assembly=([^,>]+)")


@frozen(slots=True)
class GnomadChromosomeAlleleFrequencyTask(Task):
    meta: ReferenceFileMeta
    release: GnomadRelease
    chrom: int
    fasta_task: Task
    batch_rows: int = DEFAULT_BATCH_ROWS
    max_attempts: int = DEFAULT_MAX_ATTEMPTS

    @property
    def deps(self) -> list[Task]:
        return [self.fasta_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        url = gnomad_vcf_url(self.release, self.chrom)
        fasta = load_fasta(fetch, self.fasta_task)
        assert_gnomad_header(
            read_vcf_header(url, self.max_attempts), self.release, self.chrom
        )
        tsv_path = scratch_dir / "records.tsv"
        execute_command_with_retries(
            gnomad_query_command(url, tsv_path, self.release),
            max_attempts=self.max_attempts,
        )
        out_path = scratch_dir / "allele_frequencies.parquet"
        summary = write_checked_panel(
            gnomad_tsv_batches(tsv_path, self.release, self.batch_rows),
            out_path,
            PanelBatchContext(
                fasta=fasta,
                contig_codes={gnomad_contig_name(self.release, self.chrom): self.chrom},
            ),
            byte_stream_split_columns=[
                panel_an_col(group) for group in self.release.groups
            ],
        )
        assert_gnomad_summary(summary, self.chrom)
        logger.info(
            "gnomAD chromosome allele frequencies written",
            release=self.release.name,
            chromosome=self.chrom,
            rows_written=summary.rows_written,
            fasta_ambiguous_dropped=summary.counts.fasta_ambiguous,
        )
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, release: GnomadRelease, chrom: int, fasta_task: Task
    ) -> "GnomadChromosomeAlleleFrequencyTask":
        fasta_meta = fasta_task.meta
        assert isinstance(fasta_meta, FASTAMeta), (
            f"fasta_task must carry FASTAMeta, got {type(fasta_meta).__name__}"
        )
        assert fasta_meta.build == release.build, (
            f"{release.name} is build {release.build}; the FASTA is build {fasta_meta.build}"
        )
        assert chrom in release.chromosomes, (
            f"chromosome {chrom} is not in {release.name}: {release.chromosomes}"
        )
        contig = gwaslab_code_to_contig_name(chrom)
        return cls(
            meta=ReferenceFileMeta(
                group=GNOMAD_GROUP,
                sub_group=release.name,
                sub_folder=PurePath(GNOMAD_PART_SUB_FOLDER),
                id=AssetId(f"gnomad_{release.name}_chr{contig}_allele_frequencies"),
                filename=f"chr{contig}_allele_frequencies",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            ),
            release=release,
            chrom=chrom,
            fasta_task=fasta_task,
        )


def read_vcf_header(url: str, max_attempts: int) -> str:
    return execute_command_with_retries(
        ["bcftools", "view", "-h", shlex.quote(url)], max_attempts=max_attempts
    )


def assert_gnomad_header(header: str, release: GnomadRelease, chrom: int) -> None:
    lines = header.splitlines()
    info_ids = {match.group(1) for line in lines if (match := _INFO_ID.match(line))}
    required = [
        column
        for group in release.groups
        for column in (panel_af_col(group), panel_an_col(group))
    ]
    missing = [column for column in required if column not in info_ids]
    assert not missing, f"{release.name} header lacks INFO fields {missing}"
    contig_lines = [line for line in lines if _CONTIG_ID.match(line)]
    assemblies = {
        match.group(1) for line in contig_lines if (match := _ASSEMBLY.search(line))
    }
    assert assemblies == {release.header_assembly}, (
        f"{release.name} header contig assemblies {sorted(assemblies)} != "
        f"{release.header_assembly!r}"
    )
    contigs = {
        match.group(1) for line in contig_lines if (match := _CONTIG_ID.match(line))
    }
    expected_contig = gnomad_contig_name(release, chrom)
    assert expected_contig in contigs, (
        f"{release.name} header does not declare contig {expected_contig}"
    )


def gnomad_query_command(url: str, tsv_path: Path, release: GnomadRelease) -> list[str]:
    """One bcftools query process: PASS records polymorphic in a main group, as a TSV."""
    polymorphic = " || ".join(
        f"{panel_af_col(group)}>0" for group in release.main_groups
    )
    fields = "".join(
        f"\\t%INFO/{panel_af_col(group)}\\t%INFO/{panel_an_col(group)}"
        for group in release.groups
    )
    return [
        "bcftools",
        "query",
        "-i",
        f"'FILTER=\"PASS\" && ({polymorphic})'",
        "-f",
        f"'%CHROM\\t%POS\\t%REF\\t%ALT{fields}\\n'",
        "-o",
        shlex.quote(str(tsv_path)),
        shlex.quote(url),
    ]


def gnomad_tsv_batches(
    tsv_path: Path, release: GnomadRelease, batch_rows: int
) -> Iterator[pl.DataFrame]:
    """Stream the headerless bcftools TSV in batches of batch_rows rows."""
    assert tsv_path.stat().st_size > 0, (
        f"{release.name}: no PASS records polymorphic in a main group in {tsv_path}"
    )
    schema: dict[str, pl.DataType] = {
        SOURCE_CONTIG_COL: pl.String(),
        GWASLAB_POS_COL: pl.Int32(),
        PANEL_REF_COL: pl.String(),
        PANEL_ALT_COL: pl.String(),
    }
    for group in release.groups:
        schema[panel_af_col(group)] = pl.Float32()
        schema[panel_an_col(group)] = pl.Int32()
    records = pl.scan_csv(
        tsv_path,
        separator="\t",
        has_header=False,
        schema=schema,
        null_values=[_TSV_NULL],
        quote_char=None,
    )
    # collect_batches is marked unstable in polars; the Task tests pin its batching.
    for batch in records.collect_batches(chunk_size=batch_rows):
        _assert_af_null_exactly_where_an_is_zero(batch, release.groups)
        yield batch


def _assert_af_null_exactly_where_an_is_zero(
    frame: pl.DataFrame, groups: tuple[GnomadGroup, ...]
) -> None:
    for group in groups:
        af = pl.col(panel_af_col(group))
        an = pl.col(panel_an_col(group))
        bad = frame.filter(an.is_null() | (af.is_null() != (an == 0)))
        assert bad.height == 0, (
            f"{group}: AF must be null exactly where AN is 0, and AN never null; "
            f"e.g. {bad.head(3).rows()}"
        )


def assert_gnomad_summary(summary: PanelWriteSummary, chrom: int) -> None:
    assert summary.counts.ref_mismatch == 0, (
        f"chromosome {chrom}: {summary.counts.ref_mismatch} gnomAD REF alleles differ "
        "from a pure A/C/G/T FASTA span"
    )
    assert summary.rows_written > 0, f"chromosome {chrom}: the table is empty"
    assert summary.chromosomes == frozenset({chrom}), (
        f"chromosome {chrom}: rows on {sorted(summary.chromosomes)}"
    )
