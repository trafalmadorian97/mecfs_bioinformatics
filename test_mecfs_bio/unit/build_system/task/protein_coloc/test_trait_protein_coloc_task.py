"""
Task-level tests of TraitProteinColocTask on a synthetic genome.

The world: two 4 kb chromosomes with a non-palindromic SNV every 10 bp. The trait has strong
signals at chr1:1000, chr2:2000, chr1:3000 and chr2:3500, plus one in the MHC (chr6, never fetched).
Protein OID1 shares the chr1:1000 and chr2:2000 signals and is served with half its rows
stored swapped, so orientation has to fix them. OID2 is served with too few rows to test.
OID3 is a weak candidate that the pair budget leaves out.
"""

from collections.abc import Sequence
from pathlib import Path, PurePath

import numpy as np
import polars as pl
import pytest

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
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.genome_reference_harmonization.reference_panel_task import (
    PANEL_AF_COL,
    PANEL_ALT_COL,
    PANEL_REF_COL,
)
from mecfs_bio.build_system.task.ppp_ldsc.gene_coords import (
    ST3_GENE_CHROM_COL,
    ST3_GENE_END_COL,
    ST3_GENE_START_COL,
    ST3_OLINK_ID_COL,
)
from mecfs_bio.build_system.task.protein_coloc.candidate_selection import (
    SelectionConfig,
)
from mecfs_bio.build_system.task.protein_coloc.coloc_abf import (
    QuantitativeTraitModel,
)
from mecfs_bio.build_system.task.protein_coloc.ppp_window_fetcher import LocusWindow
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    PppProteinHitTableTask,
)
from mecfs_bio.build_system.task.protein_coloc.trait_protein_coloc_task import (
    ColocRunConfig,
    TraitProteinColocTask,
)
from mecfs_bio.build_system.wf.base_wf import WF, make_wf
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
from mecfs_bio.constants.protein_coloc_constants import (
    FAILURE_REASON_COL,
    HIT_CHR_COL,
    HIT_GENE_COL,
    HIT_LEAD_MLOG10P_COL,
    HIT_LEAD_POS_COL,
    HIT_N_SIG_VARIANTS_COL,
    HIT_SEG_END_COL,
    HIT_SEG_START_COL,
    LOCUS_ID_COL,
    N_ALIGNED_COL,
    PP_H4_COL,
    PROTEIN_ID_COL,
    RAW_RESULTS_FILENAME,
    STATUS_COL,
    STATUS_EXCLUDED_MHC,
    STATUS_FAILED,
    STATUS_NO_CANDIDATES,
    STATUS_NOT_SELECTED_BUDGET,
    STATUS_TESTED,
    WINDOW_SCHEMA,
    WINDOWS_DIRNAME,
    pair_file_stem,
)
from test_mecfs_bio.unit.build_system.task.genome_reference_harmonization.genome_reference_fixtures import (
    EUR_PANEL_COLUMNS,
    write_fasta,
)

_CHROM_LENGTH = 4000
_STEP = 10
_SE = 0.02
_TRAIT_N = 50_000.0
_PROTEIN_N = 30_000.0
# chr2:3500 has no protein hit, so its locus has no candidates.
_SIGNALS = {(1, 1000): 12.0, (2, 2000): 12.0, (1, 3000): 12.0, (2, 3500): 12.0}
_PROTEIN_SIGNALS = {(1, 1000): 15.0, (2, 2000): 15.0}
_MHC_POS = 30_000_000
_TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}
_CONFIG = ColocRunConfig(
    half_window_bp=300,
    selection=SelectionConfig(max_pairs=3),
    min_aligned_variants=20,
)


def _sequences() -> dict[str, str]:
    rng = np.random.default_rng(7)
    return {f"chr{c}": "".join(rng.choice(list("ACGT"), _CHROM_LENGTH)) for c in (1, 2)}


def _variants(sequences: dict[str, str]) -> pl.DataFrame:
    """Every _STEP bp on chr1 and chr2: NEA = reference base, EA = its transition."""
    rows = [
        (
            chrom,
            pos,
            _TRANSITION[sequences[f"chr{chrom}"][pos - 1]],
            sequences[f"chr{chrom}"][pos - 1],
        )
        for chrom in (1, 2)
        for pos in range(_STEP, _CHROM_LENGTH, _STEP)
    ]
    return pl.DataFrame(
        rows,
        schema=[
            GWASLAB_CHROM_COL,
            GWASLAB_POS_COL,
            GWASLAB_EFFECT_ALLELE_COL,
            GWASLAB_NON_EFFECT_ALLELE_COL,
        ],
        orient="row",
    ).cast({GWASLAB_CHROM_COL: pl.Int64, GWASLAB_POS_COL: pl.Int64})


def _with_stats(
    variants: pl.DataFrame, signals: dict[tuple[int, int], float], n: float, seed: int
) -> pl.DataFrame:
    rng = np.random.default_rng(seed)
    z = rng.normal(0.0, 0.5, variants.height)
    for i, (chrom, pos) in enumerate(
        zip(variants[GWASLAB_CHROM_COL], variants[GWASLAB_POS_COL], strict=True)
    ):
        z[i] += signals.get((chrom, pos), 0.0)
    return variants.with_columns(
        pl.Series(GWASLAB_BETA_COL, z * _SE),
        pl.lit(_SE).alias(GWASLAB_SE_COL),
        pl.Series(
            GWASLAB_EFFECT_ALLELE_FREQ_COL, rng.uniform(0.1, 0.4, variants.height)
        ),
        pl.lit(n).alias(GWASLAB_SAMPLE_SIZE_COLUMN),
    ).cast(WINDOW_SCHEMA)


def _swap_every_other_row(frame: pl.DataFrame) -> pl.DataFrame:
    odd = pl.int_range(pl.len()) % 2 == 1
    return frame.with_columns(
        pl.when(odd)
        .then(pl.col(GWASLAB_NON_EFFECT_ALLELE_COL))
        .otherwise(pl.col(GWASLAB_EFFECT_ALLELE_COL))
        .alias(GWASLAB_EFFECT_ALLELE_COL),
        pl.when(odd)
        .then(pl.col(GWASLAB_EFFECT_ALLELE_COL))
        .otherwise(pl.col(GWASLAB_NON_EFFECT_ALLELE_COL))
        .alias(GWASLAB_NON_EFFECT_ALLELE_COL),
        pl.when(odd)
        .then(-pl.col(GWASLAB_BETA_COL))
        .otherwise(pl.col(GWASLAB_BETA_COL))
        .alias(GWASLAB_BETA_COL),
        pl.when(odd)
        .then(1 - pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
        .otherwise(pl.col(GWASLAB_EFFECT_ALLELE_FREQ_COL))
        .alias(GWASLAB_EFFECT_ALLELE_FREQ_COL),
    )


class _FakeFetcher:
    """Serves prepared dense protein data and counts fetches per protein."""

    def __init__(self, dense: dict[str, pl.DataFrame]) -> None:
        self.dense = dense
        self.calls: dict[str, int] = {}

    def fetch(
        self,
        protein_id: str,
        windows: Sequence[LocusWindow],
        scratch_dir: Path,
        wf: WF,
    ) -> dict[str, pl.DataFrame]:
        self.calls[protein_id] = self.calls.get(protein_id, 0) + 1
        frame = self.dense[protein_id]
        return {
            w.locus_id: frame.filter(
                pl.col(GWASLAB_CHROM_COL) == w.chrom,
                pl.col(GWASLAB_POS_COL).is_between(w.start, w.end),
            )
            for w in windows
        }


def _hit_segments() -> pl.DataFrame:
    rows = [
        ("OID1", 1, 1000, 100.0),
        ("OID1", 2, 2000, 100.0),
        ("OID2", 1, 3000, 80.0),
        ("OID3", 1, 1000, 9.0),
    ]
    return pl.DataFrame(
        {
            PROTEIN_ID_COL: [r[0] for r in rows],
            HIT_GENE_COL: [f"GENE_{r[0]}" for r in rows],
            HIT_CHR_COL: [r[1] for r in rows],
            HIT_SEG_START_COL: [r[2] for r in rows],
            HIT_SEG_END_COL: [r[2] for r in rows],
            HIT_LEAD_POS_COL: [r[2] for r in rows],
            HIT_LEAD_MLOG10P_COL: [r[3] for r in rows],
            HIT_N_SIG_VARIANTS_COL: [1] * len(rows),
        }
    )


def _gene_coords() -> pl.DataFrame:
    return pl.DataFrame(
        {
            ST3_OLINK_ID_COL: ["OID1"],
            ST3_GENE_CHROM_COL: ["1"],
            ST3_GENE_START_COL: ["900"],
            ST3_GENE_END_COL: ["1100"],
        }
    )


def _parquet_spec() -> DataFrameReadSpec:
    return DataFrameReadSpec(DataFrameParquetFormat())


def _table_meta(asset_id: str) -> ResultTableMeta:
    return ResultTableMeta(
        id=AssetId(asset_id),
        trait="synthetic",
        project="test",
        sub_dir=PurePath("analysis"),
        extension=".parquet",
        read_spec=_parquet_spec(),
    )


def _trait_meta(harmonized: bool) -> FilteredGWASDataMeta:
    return FilteredGWASDataMeta(
        id=AssetId("trait"),
        trait="synthetic_trait",
        project="synthetic_project",
        sub_dir="processed",
        read_spec=_parquet_spec(),
        harmonization_info=(
            HarmonizationInfo(
                build="38",
                ref_allele_col=GWASLAB_NON_EFFECT_ALLELE_COL,
                pos_col=GWASLAB_POS_COL,
            )
            if harmonized
            else None
        ),
        genome_build="38",
    )


def _fasta_task() -> FakeTask:
    return FakeTask(
        FASTAMeta(
            group="genome_sequence",
            sub_group="synthetic",
            sub_folder=PurePath("processed"),
            id=AssetId("fasta"),
            build="38",
        )
    )


def _panel_task() -> FakeTask:
    return FakeTask(
        HarmonizableReferenceTableMeta(
            group="reference_panel_allele_frequencies",
            sub_group="synthetic",
            sub_folder=PurePath("processed"),
            extension=".parquet",
            id=AssetId("panel"),
            read_spec=_parquet_spec(),
            harmonization_info=HarmonizationInfo(
                build="38", ref_allele_col=PANEL_REF_COL, pos_col=GWASLAB_POS_COL
            ),
            allele_frequency_columns=EUR_PANEL_COLUMNS,
        )
    )


def _hit_table_task() -> PppProteinHitTableTask:
    return PppProteinHitTableTask.create(
        asset_id="hits", protein_tasks=(), index_task=FakeTask(_table_meta("index"))
    )


def _create(fetcher: _FakeFetcher, trait_task: FakeTask) -> TraitProteinColocTask:
    return TraitProteinColocTask.create(
        asset_id="coloc",
        trait_task=trait_task,
        trait_model=QuantitativeTraitModel(sd_y=None),
        hit_table_task=_hit_table_task(),
        gene_coords_task=FakeTask(_table_meta("gene_coords")),
        fasta_task=_fasta_task(),
        panel_task=_panel_task(),
        panel_ancestry="eur",
        fetcher=fetcher,
        config=_CONFIG,
    )


def _write(path: Path, frame: pl.DataFrame) -> FileAsset:
    frame.write_parquet(path)
    return FileAsset(path)


def _run(tmp_path: Path) -> tuple[_FakeFetcher, Path]:
    sequences = _sequences()
    variants = _variants(sequences)
    trait = pl.concat(
        [
            _with_stats(variants, _SIGNALS, _TRAIT_N, seed=1),
            pl.DataFrame(
                {
                    GWASLAB_CHROM_COL: [6],
                    GWASLAB_POS_COL: [_MHC_POS],
                    GWASLAB_EFFECT_ALLELE_COL: ["A"],
                    GWASLAB_NON_EFFECT_ALLELE_COL: ["G"],
                    GWASLAB_BETA_COL: [20 * _SE],
                    GWASLAB_SE_COL: [_SE],
                    GWASLAB_EFFECT_ALLELE_FREQ_COL: [0.3],
                    GWASLAB_SAMPLE_SIZE_COLUMN: [_TRAIT_N],
                },
                schema=WINDOW_SCHEMA,
            ),
        ]
    )
    oid1 = _swap_every_other_row(
        _with_stats(variants, _PROTEIN_SIGNALS, _PROTEIN_N, seed=2)
    )
    oid2 = _with_stats(variants, {}, _PROTEIN_N, seed=3).filter(
        pl.col(GWASLAB_POS_COL).is_between(2990, 3030)
    )
    oid3 = _with_stats(variants, {}, _PROTEIN_N, seed=4)
    fetcher = _FakeFetcher({"OID1": oid1, "OID2": oid2, "OID3": oid3})

    fasta_dir = tmp_path / "fasta"
    write_fasta(fasta_dir, sequences)
    panel = pl.DataFrame(
        schema={
            GWASLAB_CHROM_COL: pl.Int32,
            GWASLAB_POS_COL: pl.Int32,
            PANEL_REF_COL: pl.String,
            PANEL_ALT_COL: pl.String,
            PANEL_AF_COL: pl.Float32,
        }
    )
    assets: dict[AssetId, Asset] = {
        AssetId("trait"): _write(tmp_path / "trait.parquet", trait),
        AssetId("hits"): _write(tmp_path / "hits.parquet", _hit_segments()),
        AssetId("gene_coords"): _write(
            tmp_path / "gene_coords.parquet", _gene_coords()
        ),
        AssetId("fasta"): DirectoryAsset(fasta_dir),
        AssetId("panel"): _write(tmp_path / "panel.parquet", panel),
    }
    task = _create(fetcher, FakeTask(_trait_meta(harmonized=True)))
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(
        scratch, fetch=lambda asset_id: assets[asset_id], wf=make_wf()
    )
    assert isinstance(result, DirectoryAsset)
    return fetcher, result.path


def _results(out_dir: Path) -> pl.DataFrame:
    return pl.read_parquet(out_dir / RAW_RESULTS_FILENAME)


def test_shared_signal_is_tested_with_high_h4(tmp_path: Path) -> None:
    _, out_dir = _run(tmp_path)
    rows = _results(out_dir).filter(
        pl.col(STATUS_COL) == STATUS_TESTED, pl.col(PROTEIN_ID_COL) == "OID1"
    )
    assert (rows[PP_H4_COL] > 0.8).all()
    # Orientation un-swapped the protein rows, so every window variant aligned.
    assert (rows[N_ALIGNED_COL] == 61).all()


def test_protein_at_two_loci_on_different_chromosomes_is_fetched_once(
    tmp_path: Path,
) -> None:
    fetcher, out_dir = _run(tmp_path)
    assert fetcher.calls == {"OID1": 1, "OID2": 1}
    tested = _results(out_dir).filter(
        pl.col(PROTEIN_ID_COL) == "OID1", pl.col(STATUS_COL) == STATUS_TESTED
    )
    assert tested.height == 2


def test_too_few_aligned_variants_marks_failed(tmp_path: Path) -> None:
    _, out_dir = _run(tmp_path)
    row = _results(out_dir).filter(pl.col(PROTEIN_ID_COL) == "OID2")
    assert row[STATUS_COL].item() == STATUS_FAILED
    assert row[FAILURE_REASON_COL].item() is not None
    assert row[PP_H4_COL].item() is None


def test_unselected_and_mhc_rows_are_reported(tmp_path: Path) -> None:
    _, out_dir = _run(tmp_path)
    results = _results(out_dir)
    assert (
        results.filter(pl.col(PROTEIN_ID_COL) == "OID3")[STATUS_COL].item()
        == STATUS_NOT_SELECTED_BUDGET
    )
    mhc = results.filter(pl.col(STATUS_COL) == STATUS_EXCLUDED_MHC)
    assert mhc.height == 1
    assert mhc[PROTEIN_ID_COL].item() is None


def test_integer_columns_are_int32(tmp_path: Path) -> None:
    _, out_dir = _run(tmp_path)
    schema = _results(out_dir).schema
    assert pl.Int64 not in schema.dtypes()
    assert schema[N_ALIGNED_COL] == pl.Int32


def test_window_parquet_written_for_each_tested_pair(tmp_path: Path) -> None:
    _, out_dir = _run(tmp_path)
    tested = _results(out_dir).filter(pl.col(STATUS_COL) == STATUS_TESTED)
    for locus_id, protein_id in tested.select(LOCUS_ID_COL, PROTEIN_ID_COL).rows():
        path = (
            out_dir
            / WINDOWS_DIRNAME
            / f"{pair_file_stem(locus_id, protein_id)}.parquet"
        )
        assert pl.read_parquet(path).height == 61


def test_create_rejects_unharmonized_trait() -> None:
    with pytest.raises(AssertionError):
        _create(_FakeFetcher({}), FakeTask(_trait_meta(harmonized=False)))


def test_locus_without_candidates_is_reported(tmp_path: Path) -> None:
    _, out_dir = _run(tmp_path)
    rows = _results(out_dir).filter(pl.col(STATUS_COL) == STATUS_NO_CANDIDATES)
    assert rows[LOCUS_ID_COL].to_list() == ["chr2:3200-3800"]
    assert rows[PROTEIN_ID_COL].item() is None
