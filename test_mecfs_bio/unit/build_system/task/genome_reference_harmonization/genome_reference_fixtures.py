"""Synthetic genome, reference panel and Task runner for genome-reference harmonization tests."""

from collections.abc import Mapping, Sequence
from pathlib import Path, PurePath

import polars as pl
import pysam
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.fasta_meta import FASTAMeta
from mecfs_bio.build_system.meta.reference_meta.harmonizable_reference_table_meta import (
    HarmonizableReferenceTableMeta,
)
from mecfs_bio.build_system.meta.reference_meta.panel_allele_frequency_columns import (
    PanelAlleleFrequencyColumns,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    FASTA_FILENAME,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.allele_frequency_panel_constants import (
    PanelAncestry,
    panel_af_col,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
)

# Line width 7 so that alleles cross FASTA line breaks. 1-based landmarks:
#   1-10  A C G T T G G G G G   (5: T followed by a G run -> T/TG is ambiguous)
#   11-20 C A T G C A A T T C   (13, 18, 19: T for A/T palindromes)
#   21-30 G A T C G A T C G A
#   31-40 T x 10                (a homopolymer for T/TT ambiguous indels)
#   41-50 A C A C A C A C A C   (41: AC for the AC/GT palindromic MNP)
#   51-60 G T C A G T C A G T
LINE_WIDTH = 7
CHR1_SEQUENCE = "ACGTTGGGGGCATGCAATTCGATCGATCGATTTTTTTTTTACACACACACGTCAGTCAGT"
CHR2_SEQUENCE = "GGGGAAAACCCCTTTT"
SE_VALUE = 0.01

# min_checkable_ambiguous_indels=1 so a single suspicious ambiguous indel can drive the
# suspicious-fraction trust gate in tests.
TEST_OPTIONS = GenomeReferenceHarmonizationOptions(
    min_checkable_snvs=1, min_checkable_indels=1, min_checkable_ambiguous_indels=1
)

_SUMSTATS_ID = "sumstats"
_FASTA_ID = "fasta"
_PANEL_ID = "panel"


@frozen(slots=True)
class Variant:
    pos: int
    ea: str
    nea: str
    eaf: float | None = 0.3
    beta: float = 0.1
    chrom: int = 1


@frozen(slots=True)
class PanelRecord:
    pos: int
    ref: str
    alt: str
    af: float
    chrom: int = 1


# A table of only these two is fully reference-consistent (trusted under TEST_OPTIONS);
# adding INCONSISTENT_SNV (EA is the reference base) makes it untrusted.
CONSISTENT_SNV = Variant(pos=1, ea="G", nea="A")
CONSISTENT_INDEL = Variant(pos=21, ea="GC", nea="G")
INCONSISTENT_SNV = Variant(pos=2, ea="C", nea="T", beta=0.2)


def sumstats_frame(variants: Sequence[Variant]) -> pl.DataFrame:
    return pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [v.chrom for v in variants],
            GWASLAB_POS_COL: [v.pos for v in variants],
            GWASLAB_EFFECT_ALLELE_COL: [v.ea for v in variants],
            GWASLAB_NON_EFFECT_ALLELE_COL: [v.nea for v in variants],
            GWASLAB_EFFECT_ALLELE_FREQ_COL: [v.eaf for v in variants],
            GWASLAB_BETA_COL: [v.beta for v in variants],
            GWASLAB_SE_COL: [SE_VALUE for _ in variants],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int64,
            GWASLAB_POS_COL: pl.Int64,
            GWASLAB_EFFECT_ALLELE_COL: pl.String,
            GWASLAB_NON_EFFECT_ALLELE_COL: pl.String,
            GWASLAB_EFFECT_ALLELE_FREQ_COL: pl.Float64,
            GWASLAB_BETA_COL: pl.Float64,
            GWASLAB_SE_COL: pl.Float64,
        },
    )


SYNTHETIC_SEQUENCES: Mapping[str, str] = {"chr1": CHR1_SEQUENCE, "chr2": CHR2_SEQUENCE}


def write_fasta(
    directory: Path, sequences: Mapping[str, str] = SYNTHETIC_SEQUENCES
) -> None:
    """Write an uncompressed, faidx-indexed FASTA with LINE_WIDTH-base lines."""
    directory.mkdir()
    fasta_path = directory / FASTA_FILENAME
    lines: list[str] = []
    for name, sequence in sequences.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[start : start + LINE_WIDTH]
            for start in range(0, len(sequence), LINE_WIDTH)
        )
    fasta_path.write_text("\n".join(lines) + "\n")
    pysam.faidx(str(fasta_path))


def _write_panel(path: Path, records: Sequence[PanelRecord]) -> None:
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [r.chrom for r in records],
            GWASLAB_POS_COL: [r.pos for r in records],
            PANEL_REF_COL: [r.ref for r in records],
            PANEL_ALT_COL: [r.alt for r in records],
            PANEL_AF_COL: [r.af for r in records],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            PANEL_AF_COL: pl.Float32,
        },
    ).sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL).write_parquet(path)


EUR_PANEL_COLUMNS = PanelAlleleFrequencyColumns.single("eur", PANEL_AF_COL)
ANCESTRY_PANEL_COLUMNS = PanelAlleleFrequencyColumns.prefixed(("nfe", "nfe_nwe"))


@frozen(slots=True)
class AncestryPanelRecord:
    """A panel record with gnomAD-style AF_nfe and AF_nfe_nwe columns (None = AN 0)."""

    pos: int
    ref: str
    alt: str
    nfe: float | None
    nfe_nwe: float | None
    chrom: int = 1


def _write_ancestry_panel(path: Path, records: Sequence[AncestryPanelRecord]) -> None:
    pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [r.chrom for r in records],
            GWASLAB_POS_COL: [r.pos for r in records],
            PANEL_REF_COL: [r.ref for r in records],
            PANEL_ALT_COL: [r.alt for r in records],
            panel_af_col("nfe"): [r.nfe for r in records],
            panel_af_col("nfe_nwe"): [r.nfe_nwe for r in records],
        },
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            panel_af_col("nfe"): pl.Float32,
            panel_af_col("nfe_nwe"): pl.Float32,
        },
    ).sort(GWASLAB_CHROM_COL, GWASLAB_POS_COL).write_parquet(path)


def harmonization_task(
    panel_columns: PanelAlleleFrequencyColumns | None,
    panel_ancestry: PanelAncestry,
    options: GenomeReferenceHarmonizationOptions = TEST_OPTIONS,
    pipe: DataProcessingPipe = IdentityPipe(),
) -> GenomeReferenceHarmonizationTask:
    """The Task under test, wired to FakeTask inputs whose ids run_harmonization serves."""
    parquet_spec = DataFrameReadSpec(DataFrameParquetFormat())
    return GenomeReferenceHarmonizationTask.create(
        asset_id="harmonized",
        sumstats_task=FakeTask(
            FilteredGWASDataMeta(
                id=AssetId(_SUMSTATS_ID),
                trait="synthetic_trait",
                project="synthetic_project",
                sub_dir="processed",
                read_spec=parquet_spec,
            )
        ),
        fasta_task=FakeTask(
            FASTAMeta(
                group="genome_sequence",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                id=AssetId(_FASTA_ID),
                build="19",
            )
        ),
        panel_task=FakeTask(
            HarmonizableReferenceTableMeta(
                group="reference_panel_allele_frequencies",
                sub_group="synthetic",
                sub_folder=PurePath("processed"),
                extension=".parquet",
                id=AssetId(_PANEL_ID),
                read_spec=parquet_spec,
                harmonization_info=HarmonizationInfo(
                    build="19", ref_allele_col=PANEL_REF_COL, pos_col=GWASLAB_POS_COL
                ),
                allele_frequency_columns=panel_columns,
            )
        ),
        panel_ancestry=panel_ancestry,
        options=options,
        pipe=pipe,
    )


def run_harmonization(
    work_dir: Path,
    sumstats: pl.DataFrame,
    panel: Sequence[PanelRecord] = (),
    options: GenomeReferenceHarmonizationOptions = TEST_OPTIONS,
    pipe: DataProcessingPipe = IdentityPipe(),
    ancestry_panel: Sequence[AncestryPanelRecord] | None = None,
    panel_ancestry: PanelAncestry = "eur",
) -> pl.DataFrame:
    """Execute the Task on synthetic inputs in a fresh work_dir and return the output table.

    By default the panel is the single-ancestry "eur" panel built from panel. Pass
    ancestry_panel to use a panel with AF_nfe and AF_nfe_nwe columns instead.
    """
    work_dir.mkdir(parents=True)
    sumstats_path = work_dir / "sumstats.parquet"
    sumstats.write_parquet(sumstats_path)
    fasta_dir = work_dir / "fasta"
    write_fasta(fasta_dir)
    panel_path = work_dir / "panel.parquet"
    if ancestry_panel is None:
        _write_panel(panel_path, panel)
        panel_columns = EUR_PANEL_COLUMNS
    else:
        _write_ancestry_panel(panel_path, ancestry_panel)
        panel_columns = ANCESTRY_PANEL_COLUMNS
    task = harmonization_task(
        panel_columns=panel_columns,
        panel_ancestry=panel_ancestry,
        options=options,
        pipe=pipe,
    )
    assets: dict[str, Asset] = {
        _SUMSTATS_ID: FileAsset(sumstats_path),
        _FASTA_ID: DirectoryAsset(fasta_dir),
        _PANEL_ID: FileAsset(panel_path),
    }

    def fetch(asset_id: AssetId) -> Asset:
        return assets[asset_id]

    scratch = work_dir / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, FileAsset)
    return pl.read_parquet(result.path)


def positions(frame: pl.DataFrame, chrom: int = 1) -> list[int]:
    return frame.filter(pl.col(GWASLAB_CHROM_COL) == chrom)[GWASLAB_POS_COL].to_list()


def row_at(frame: pl.DataFrame, pos: int, chrom: int = 1) -> dict[str, object]:
    matching = frame.filter(
        (pl.col(GWASLAB_CHROM_COL) == chrom) & (pl.col(GWASLAB_POS_COL) == pos)
    )
    assert matching.height == 1, (
        f"expected one row at {chrom}:{pos}, got {matching.height}"
    )
    return matching.row(0, named=True)
