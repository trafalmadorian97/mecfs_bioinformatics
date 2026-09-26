"""
Convert a reference panel VCF into a sites-only allele-frequency parquet.

Genome-reference harmonization looks up panel allele frequencies for palindromic SNVs
and for indels whose alleles both match the genome. This Task strips genotypes with
bcftools and keeps CHR (gwaslab numeric coding), POS, REF, ALT and AF.

- Only main contigs are kept (1-22, X, Y, MT).
- Records with no ALT or no AF are dropped.
- The panel is expected to be split and normalized already. A multi-allelic record's
  comma-separated AF fails the Float32 parse, which stops the build.
- Exact duplicate (CHR, POS, REF, ALT) keys with one AF collapse to one row. Keys with
  conflicting AF are removed and counted in the log, so a lookup never sees two answers.

Rows are sorted by CHR, POS so a per-chromosome scan reads only matching row groups.
The work is done one contig at a time to bound memory. Wrap instances in
DiscardDepsWrapper so the multi-gigabyte VCF is not also kept in the asset store.
bcftools is resolved from the pixi environment.
"""

from pathlib import Path, PurePath

import polars as pl
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    contig_to_gwaslab_code,
)
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL
from mecfs_bio.util.subproc.run_command import execute_command

logger = structlog.get_logger()

PANEL_REF_COL = "REF"
PANEL_ALT_COL = "ALT"
PANEL_AF_COL = "AF"
PANEL_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    PANEL_REF_COL,
    PANEL_ALT_COL,
    PANEL_AF_COL,
]
_CONTIG_COL = "contig"
_N_AF_COL = "n_distinct_af"
_SITE_KEYS = [GWASLAB_POS_COL, PANEL_REF_COL, PANEL_ALT_COL]


def _write_sites_tsv(vcf_path: Path, tsv_path: Path) -> None:
    execute_command(
        [
            "bcftools",
            "view",
            "-G",
            "-Ou",
            str(vcf_path),
            "|",
            "bcftools",
            "query",
            "-f",
            r"'%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\n'",
            "-o",
            str(tsv_path),
        ]
    )


def _scan_sites(tsv_path: Path) -> pl.LazyFrame:
    return pl.scan_csv(
        tsv_path,
        separator="\t",
        has_header=False,
        null_values=["."],
        schema={
            _CONTIG_COL: pl.String,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            PANEL_AF_COL: pl.Float32,
        },
    )


def _contigs_by_code(sites: pl.LazyFrame) -> dict[int, list[str]]:
    names = sites.select(pl.col(_CONTIG_COL).unique()).collect(engine="streaming")[
        _CONTIG_COL
    ]
    by_code: dict[int, list[str]] = {}
    skipped: list[str] = []
    for name in names.to_list():
        code = contig_to_gwaslab_code(name)
        if code is None:
            skipped.append(name)
        else:
            by_code.setdefault(code, []).append(name)
    if skipped:
        logger.info("skipping non-main panel contigs", contigs=sorted(skipped))
    return by_code


def _one_chromosome(sites: pl.LazyFrame, code: int, names: list[str]) -> pl.DataFrame:
    rows = (
        sites.filter(pl.col(_CONTIG_COL).is_in(names))
        .drop_nulls([PANEL_ALT_COL, PANEL_AF_COL])
        .collect(engine="streaming")
    )
    grouped = rows.group_by(_SITE_KEYS).agg(
        pl.col(PANEL_AF_COL).first(),
        pl.col(PANEL_AF_COL).n_unique().alias(_N_AF_COL),
    )
    n_conflicting = grouped.filter(pl.col(_N_AF_COL) > 1).height
    if n_conflicting:
        logger.warning(
            "removing panel sites with conflicting duplicate allele frequencies",
            chromosome=code,
            n_sites=n_conflicting,
        )
    return (
        grouped.filter(pl.col(_N_AF_COL) == 1)
        .with_columns(pl.lit(code, dtype=pl.Int32).alias(GWASLAB_CHROM_COL))
        .select(PANEL_COLUMNS)
        .sort(_SITE_KEYS)
    )


@frozen(slots=True)
class ReferencePanelAlleleFrequencyTask(Task):
    meta: HarmonizableReferenceTableMeta
    vcf_task: Task

    @property
    def deps(self) -> list[Task]:
        return [self.vcf_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        vcf = fetch(self.vcf_task.asset_id)
        assert isinstance(vcf, FileAsset), (
            f"expected {self.vcf_task.asset_id} to be a FileAsset, got {type(vcf).__name__}"
        )
        tsv_path = scratch_dir / "panel_sites.tsv"
        _write_sites_tsv(vcf_path=vcf.path, tsv_path=tsv_path)
        sites = _scan_sites(tsv_path)
        parts_dir = scratch_dir / "parts"
        parts_dir.mkdir()
        part_paths: list[Path] = []
        for code, names in sorted(_contigs_by_code(sites).items()):
            part_path = parts_dir / f"chr{code}.parquet"
            _one_chromosome(sites, code, names).write_parquet(part_path)
            part_paths.append(part_path)
        assert part_paths, f"no main-contig records in {vcf.path}"
        out_path = scratch_dir / "panel_allele_frequencies.parquet"
        pl.concat([pl.scan_parquet(path) for path in part_paths]).sink_parquet(out_path)
        return FileAsset(out_path)

    @classmethod
    def create(
        cls, vcf_task: Task, asset_id: str, build: GenomeBuild
    ) -> "ReferencePanelAlleleFrequencyTask":
        source_meta = vcf_task.meta
        assert isinstance(source_meta, ReferenceFileMeta), (
            f"expected a ReferenceFileMeta source for {asset_id}, got {type(source_meta).__name__}"
        )
        return cls(
            meta=HarmonizableReferenceTableMeta(
                group="reference_panel_allele_frequencies",
                sub_group=source_meta.sub_group,
                sub_folder=PurePath("processed"),
                id=AssetId(asset_id),
                filename="panel_allele_frequencies",
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
                harmonization_info=HarmonizationInfo(
                    build=build, ref_allele_col=PANEL_REF_COL, pos_col=GWASLAB_POS_COL
                ),
            ),
            vcf_task=vcf_task,
        )
