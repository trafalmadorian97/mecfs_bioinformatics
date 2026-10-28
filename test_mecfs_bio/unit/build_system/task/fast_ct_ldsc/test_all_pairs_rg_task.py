"""The all-pairs rg and heritability tasks on synthetic traits with known genetic correlation.

Four traits on S = 20_000 SNPs, simulated under the LDSC model z_t = sqrt(N_t h2_t l / M) g_t + e_t
with independent noise e (no sample overlap): A and B share genetic effects with rg = 0.6, C is
genetically independent of both, and D has chi-square falling with LD score, so its estimated h2
is negative. B carries a per-SNP N; A is binary.
"""

from pathlib import Path

import numpy as np
import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fast_ct_ldsc.aligned_trait_task import (
    FastCtLdscAlignedTraitTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import AlignmentGuards
from mecfs_bio.build_system.task.fast_ct_ldsc.all_pairs_rg_task import (
    FastCtLdscAllPairsRgTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.heritability_task import (
    FastCtLdscHeritabilityTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.trait_matrix import load_trait_matrix
from mecfs_bio.build_system.task.genetic_correlation_clustermap_task import (
    XR_GENETIC_CORR_ARRAY,
    GeneticCorrSource,
    load_xr_corr_dataset,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    PhenotypeInfo,
    QuantPhenotype,
    SumstatsSource,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import liability_h2
from mecfs_bio.build_system.task.lava_task import _get_sample_overlap_path
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    FAST_CT_LDSC_ALIGNED_FILENAME,
    FAST_CT_LDSC_DIAGNOSTICS_FILENAME,
    FAST_CT_LDSC_H2_INT_COL,
    FAST_CT_LDSC_H2_LIAB_COL,
    FAST_CT_LDSC_H2_LIAB_SE_COL,
    FAST_CT_LDSC_H2_OBS_COL,
    FAST_CT_LDSC_H2_OBS_P1_COL,
    FAST_CT_LDSC_H2_OBS_P2_COL,
    FAST_CT_LDSC_H2_OBS_SE_COL,
    FAST_CT_LDSC_INTERCEPT_COL,
    FAST_CT_LDSC_P1_COL,
    FAST_CT_LDSC_P2_COL,
    FAST_CT_LDSC_RG_COL,
    FAST_CT_LDSC_RG_P_COL,
    FAST_CT_LDSC_RG_SE_COL,
    FAST_CT_LDSC_TRAIT_COL,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL
from test_mecfs_bio.unit.build_system.task.fast_ct_ldsc.fast_ct_ldsc_fixtures import (
    SyntheticIndex,
    build_synthetic_index,
    parquet_trait_meta,
    register_frame,
    run_task,
    trait_frame,
)

_S = 20_000
_N = 10_000.0
_H2 = 0.3
_RG_AB = 0.6
_N_BLOCKS = 200
_INFO_A = BinaryPhenotypeSampleInfo(
    sample_prevalence=0.3, estimated_population_prevalence=0.05
)


@frozen(slots=True)
class _Run:
    idx: SyntheticIndex
    aligned: dict[str, FastCtLdscAlignedTraitTask]
    pairs_task: FastCtLdscAllPairsRgTask
    heritability_task: FastCtLdscHeritabilityTask
    pairs: pl.DataFrame
    heritability: pl.DataFrame
    tmp_path: Path

    def fetch(self, asset_id: AssetId) -> Asset:
        return self.idx.assets[asset_id]


def _simulate(ld: np.ndarray) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Index-oriented z and per-SNP N for each trait."""
    rng = np.random.default_rng(11)
    s = ld.shape[0]
    m = float(s)
    g_a = rng.standard_normal(s)
    g_b = _RG_AB * g_a + np.sqrt(1.0 - _RG_AB**2) * rng.standard_normal(s)
    g_c = rng.standard_normal(s)
    n_b = rng.uniform(8_000.0, 12_000.0, size=s)

    def z(*, g: np.ndarray, n: np.ndarray | float) -> np.ndarray:
        return np.sqrt(n * _H2 * ld / m) * g + rng.standard_normal(s)

    n_const = np.full(s, _N)
    z_d = rng.standard_normal(s) * np.sqrt(np.maximum(0.2, 1.5 - 0.02 * ld))
    return {
        "A": (z(g=g_a, n=_N), n_const),
        "B": (z(g=g_b, n=n_b), n_b),
        "C": (z(g=g_c, n=_N), n_const),
        "D": (z_d, n_const),
    }


def _sample_info(alias: str) -> PhenotypeInfo:
    return _INFO_A if alias == "A" else QuantPhenotype()


def _build(tmp_path: Path) -> _Run:
    idx = build_synthetic_index(tmp_path, _S, seed=3)
    ld = idx.variants["L2"].to_numpy()
    aligned: dict[str, FastCtLdscAlignedTraitTask] = {}
    for alias, (z, n) in _simulate(ld).items():
        trait_task = FakeTask(
            meta=parquet_trait_meta(asset_id=f"trait_{alias}", genome_build="19")
        )
        register_frame(
            tmp_path, idx.assets, trait_task, trait_frame(idx.variants, z, n=n)
        )
        source = SumstatsSource(
            task=trait_task, alias=alias, sample_info=_sample_info(alias)
        )
        task = FastCtLdscAlignedTraitTask.create(
            f"aligned_{alias}",
            source=source,
            index_task=idx.index_task,
            guards=AlignmentGuards(),
        )
        idx.assets[task.asset_id] = run_task(tmp_path, task, idx.assets)
        aligned[alias] = task
    tasks = tuple(aligned.values())
    pairs_task = FastCtLdscAllPairsRgTask.create(
        "pairs",
        index_task=idx.index_task,
        aligned_tasks=tasks,
        n_blocks=_N_BLOCKS,
        batch_size=2,
    )
    heritability_task = FastCtLdscHeritabilityTask.create(
        "heritability",
        index_task=idx.index_task,
        aligned_tasks=tasks,
        n_blocks=_N_BLOCKS,
    )
    tables = {}
    for task in (pairs_task, heritability_task):
        asset = run_task(tmp_path, task, idx.assets)
        idx.assets[task.asset_id] = asset
        assert isinstance(asset, FileAsset)
        tables[task.asset_id] = pl.read_parquet(asset.path)
    return _Run(
        idx=idx,
        aligned=aligned,
        pairs_task=pairs_task,
        heritability_task=heritability_task,
        pairs=tables[pairs_task.asset_id],
        heritability=tables[heritability_task.asset_id],
        tmp_path=tmp_path,
    )


@pytest.fixture(scope="module")
def run(tmp_path_factory: pytest.TempPathFactory) -> _Run:
    return _build(tmp_path_factory.mktemp("all_pairs"))


def _row(df: pl.DataFrame, *, p1: str, p2: str) -> dict:
    rows = df.filter(
        (pl.col(FAST_CT_LDSC_P1_COL) == p1) & (pl.col(FAST_CT_LDSC_P2_COL) == p2)
    ).to_dicts()
    assert len(rows) == 1
    return rows[0]


def _herit(run: _Run, alias: str) -> dict:
    rows = run.heritability.filter(pl.col(FAST_CT_LDSC_TRAIT_COL) == alias).to_dicts()
    assert len(rows) == 1
    return rows[0]


def test_all_pairs_outputs_each_unordered_pair_once(run: _Run) -> None:
    pairs = list(zip(run.pairs[FAST_CT_LDSC_P1_COL], run.pairs[FAST_CT_LDSC_P2_COL]))
    assert sorted(pairs) == [
        ("A", "B"),
        ("A", "C"),
        ("A", "D"),
        ("B", "C"),
        ("B", "D"),
        ("C", "D"),
    ]


def test_correlated_pair_recovers_rg(run: _Run) -> None:
    ab = _row(run.pairs, p1="A", p2="B")
    assert ab[FAST_CT_LDSC_RG_COL] == pytest.approx(_RG_AB, abs=0.15)
    assert ab[FAST_CT_LDSC_RG_P_COL] < 1e-6
    assert abs(_row(run.pairs, p1="A", p2="C")[FAST_CT_LDSC_RG_COL]) < 0.2
    assert _herit(run, "A")[FAST_CT_LDSC_H2_OBS_COL] == pytest.approx(_H2, abs=0.1)


def test_nonpositive_h2_trait_gives_nan_rg_without_crash(run: _Run) -> None:
    assert _herit(run, "D")[FAST_CT_LDSC_H2_OBS_COL] <= 0.0
    for p1, p2 in (("A", "D"), ("B", "D"), ("C", "D")):
        row = _row(run.pairs, p1=p1, p2=p2)
        assert np.isnan(row[FAST_CT_LDSC_RG_COL])
        assert np.isnan(row[FAST_CT_LDSC_RG_SE_COL])


def test_pairs_and_heritability_tables_agree(run: _Run) -> None:
    for row in run.pairs.to_dicts():
        h1 = _herit(run, row[FAST_CT_LDSC_P1_COL])
        h2 = _herit(run, row[FAST_CT_LDSC_P2_COL])
        assert row[FAST_CT_LDSC_H2_OBS_P1_COL] == pytest.approx(
            h1[FAST_CT_LDSC_H2_OBS_COL], rel=1e-10
        )
        assert row[FAST_CT_LDSC_H2_OBS_P2_COL] == pytest.approx(
            h2[FAST_CT_LDSC_H2_OBS_COL], rel=1e-10
        )
        assert row[FAST_CT_LDSC_H2_INT_COL] == pytest.approx(
            h2[FAST_CT_LDSC_INTERCEPT_COL], rel=1e-10
        )


def test_liability_columns_only_for_binary(run: _Run) -> None:
    a = _herit(run, "A")
    factor = liability_h2(1.0, _INFO_A)
    assert factor is not None
    assert a[FAST_CT_LDSC_H2_LIAB_COL] == pytest.approx(
        factor * a[FAST_CT_LDSC_H2_OBS_COL], rel=1e-10
    )
    assert a[FAST_CT_LDSC_H2_LIAB_SE_COL] == pytest.approx(
        factor * a[FAST_CT_LDSC_H2_OBS_SE_COL], rel=1e-10
    )
    c = _herit(run, "C")
    assert c[FAST_CT_LDSC_H2_LIAB_COL] is None
    assert c[FAST_CT_LDSC_H2_LIAB_SE_COL] is None


def test_outputs_feed_clustermap_and_lava_overlap_readers(
    run: _Run, tmp_path: Path
) -> None:
    ds = load_xr_corr_dataset(GeneticCorrSource(task=run.pairs_task), run.fetch)
    assert ds[XR_GENETIC_CORR_ARRAY].shape == (4, 4)
    path = _get_sample_overlap_path(
        ct_ldsc_task=run.pairs_task,
        heritability_task=run.heritability_task,
        fetch=run.fetch,
        pheno_names=["A", "B", "C"],
        tmp_dir=tmp_path,
    )
    assert path is not None


def test_aligned_trait_from_a_different_index_rejected(
    run: _Run, tmp_path: Path
) -> None:
    good = run.aligned["A"]
    good_asset = run.fetch(good.asset_id)
    assert isinstance(good_asset, DirectoryAsset)
    bad_dir = tmp_path / "bad"
    bad_dir.mkdir()
    frame = pl.read_parquet(good_asset.path / FAST_CT_LDSC_ALIGNED_FILENAME)
    frame.with_columns(pl.col(GWASLAB_RSID_COL).reverse()).write_parquet(
        bad_dir / FAST_CT_LDSC_ALIGNED_FILENAME
    )
    (bad_dir / FAST_CT_LDSC_DIAGNOSTICS_FILENAME).write_text(
        (good_asset.path / FAST_CT_LDSC_DIAGNOSTICS_FILENAME).read_text()
    )
    bad = run.aligned["B"]
    assets = dict(run.idx.assets)
    assets[bad.asset_id] = DirectoryAsset(bad_dir)
    with pytest.raises(AssertionError):
        load_trait_matrix(
            run.idx.index_task, (good, bad), lambda asset_id: assets[asset_id]
        )


def test_duplicate_aliases_rejected(run: _Run) -> None:
    a = run.aligned["A"]
    with pytest.raises(AssertionError):
        FastCtLdscAllPairsRgTask.create(
            "dup",
            index_task=run.idx.index_task,
            aligned_tasks=(a, a),
            n_blocks=_N_BLOCKS,
            batch_size=2,
        )
