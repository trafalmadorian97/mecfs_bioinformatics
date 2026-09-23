import json
from pathlib import Path, PurePath

import numpy as np
import polars as pl
from attrs import frozen
from scipy.stats import spearmanr

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.reference_meta.reference_data_directory_meta import (
    ReferenceDataDirectoryMeta,
)
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.annotation_weights.l2_sldsc_snpvar_task import (
    DIAGNOSTICS_JSON_FILENAME,
    SNPVAR_COL,
    SNPVAR_PARQUET_FILENAME,
    TAU_EVEN_WEIGHTS_FILENAME,
    TAU_ODD_WEIGHTS_FILENAME,
    L2RegularizedSldscSnpvarTask,
)
from mecfs_bio.build_system.task.annotation_weights.ridge_annotation_weights_task import (
    ANNOTATION_COL,
    FAMILY_COL,
    GAMMA_RAW_COL,
)
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_SE_COL,
    GWASLAB_Z_COL,
)
from mecfs_bio.constants.polyfun_constants import (
    POLYFUN_A1_COL,
    POLYFUN_A2_COL,
    POLYFUN_BP_COL,
    POLYFUN_CHR_COL,
    POLYFUN_SNP_COL,
    POLYFUN_TO_GWASLAB_KEY_RENAME,
)

ENRICHED = "Coding_UCSC_common"
ANNOTS = [ENRICHED, "Repressed_Hoffman_common"]
MAFBINS = [f"MAFbin_lowfreq_{i}" for i in range(1, 11)] + [
    f"MAFbin_frequent_{i}" for i in range(1, 11)
]
N_CHROM = 6
N_PER_CHROM = 300
EFFECTIVE_N = 100_000
# Variants at these per-chromosome indices have no sumstats: they are scored but
# never enter the regression.
NO_SUMSTATS_INDICES = range(0, 10)
_KEY = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_COL,
]


@frozen
class _Fixture:
    task: L2RegularizedSldscSnpvarTask
    ldscore_dir: Path
    annot_path: Path
    sumstats_path: Path

    def fetch(self, asset_id: AssetId) -> Asset:
        if asset_id == self.task.annotation_ldscore_members_task.asset_id:
            return DirectoryAsset(self.ldscore_dir)
        if asset_id == self.task.annotation_matrix_task.asset_id:
            return FileAsset(self.annot_path)
        if asset_id == self.task.sumstats_task.asset_id:
            return FileAsset(self.sumstats_path)
        raise ValueError(f"unknown asset id {asset_id}")


def _build_fixture(
    tmp_path: Path, perturb_chrom: int | None, z_as_beta_se: bool = False
) -> _Fixture:
    """Synthetic data in which the ENRICHED annotation carries most heritability.

    Each variant's annotation LD-score and annotation value are the same number, so
    chi2_i = 1 + N * sum_c ell(i, c) tau_c. Each variant sits in exactly one MAF
    bin, so the MAFbin LD-scores sum to its total LD-score l_i. The MAFbins carry a
    small uniform tau so that chi2 rises with l, as in real data; otherwise h2bar
    would sit at its floor and the regression weights would ignore the sumstats.
    perturb_chrom, if given, rescales that chromosome's Z. z_as_beta_se replaces
    the Z column with BETA and SE columns whose ratio is Z.
    """
    tmp_path.mkdir(parents=True)
    rng = np.random.default_rng(0)
    cols = ANNOTS + MAFBINS
    tau = np.zeros(len(cols))
    tau[0] = 2e-5
    tau[len(ANNOTS) :] = 5e-7
    ldscore_dir = tmp_path / "ldscores"
    ldscore_dir.mkdir()
    annot_frames: list[pl.DataFrame] = []
    sumstats_frames: list[pl.DataFrame] = []
    for chrom in range(1, N_CHROM + 1):
        values = np.zeros((N_PER_CHROM, len(cols)))
        values[:, : len(ANNOTS)] = rng.uniform(0, 1, size=(N_PER_CHROM, len(ANNOTS)))
        which_bin = rng.integers(0, len(MAFBINS), size=N_PER_CHROM)
        values[np.arange(N_PER_CHROM), len(ANNOTS) + which_bin] = rng.uniform(
            1, 5, size=N_PER_CHROM
        )
        chi2 = (
            1.0
            + EFFECTIVE_N * (values @ tau)
            + rng.normal(scale=0.05, size=N_PER_CHROM)
        )
        z = np.sqrt(np.maximum(chi2, 1e-6))
        if chrom == perturb_chrom:
            z = z * 1.5
        key = pl.DataFrame(
            {
                POLYFUN_CHR_COL: [chrom] * N_PER_CHROM,
                POLYFUN_BP_COL: np.arange(1, N_PER_CHROM + 1),
                POLYFUN_SNP_COL: [f"rs{chrom}_{i}" for i in range(N_PER_CHROM)],
                POLYFUN_A1_COL: ["A"] * N_PER_CHROM,
                POLYFUN_A2_COL: ["G"] * N_PER_CHROM,
            }
        )
        with_values = key.with_columns(
            pl.Series(name, values[:, k]) for k, name in enumerate(cols)
        )
        with_values.write_parquet(
            ldscore_dir / f"baselineLF2.2.UKB.{chrom}.l2.ldscore.parquet"
        )
        annot_frames.append(with_values)
        keep = np.ones(N_PER_CHROM, dtype=bool)
        keep[list(NO_SUMSTATS_INDICES)] = False
        sumstats_frames.append(
            key.rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
            .drop(POLYFUN_SNP_COL)
            .with_columns(pl.Series(GWASLAB_Z_COL, z))
            .filter(pl.Series(keep))
        )
    annot_path = tmp_path / "annot.parquet"
    pl.concat(annot_frames).write_parquet(annot_path)
    sumstats_path = tmp_path / "sumstats.parquet"
    sumstats = pl.concat(sumstats_frames)
    if z_as_beta_se:
        se = pl.Series(GWASLAB_SE_COL, rng.uniform(0.01, 0.1, size=sumstats.height))
        sumstats = (
            sumstats.with_columns(se)
            .with_columns(
                (pl.col(GWASLAB_Z_COL) * pl.col(GWASLAB_SE_COL)).alias(GWASLAB_BETA_COL)
            )
            .drop(GWASLAB_Z_COL)
        )
    sumstats.write_parquet(sumstats_path)

    task = L2RegularizedSldscSnpvarTask.create(
        asset_id="snpvar",
        sumstats_task=FakeTask(
            FilteredGWASDataMeta(
                id=AssetId("sumstats"),
                trait="trait",
                project="project",
                sub_dir=PurePath("processed"),
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            )
        ),
        effective_sample_size=EFFECTIVE_N,
        annotation_ldscore_members_task=FakeTask(
            ReferenceDataDirectoryMeta(
                group="polyfun",
                sub_group="annotations",
                sub_folder=PurePath("raw"),
                id=AssetId("ldscores"),
            )
        ),
        annotation_matrix_task=FakeTask(
            ReferenceFileMeta(
                group="polyfun",
                sub_group="annotations",
                sub_folder=PurePath("raw"),
                id=AssetId("annot"),
                extension=".parquet",
                read_spec=DataFrameReadSpec(DataFrameParquetFormat()),
            )
        ),
    )
    return _Fixture(
        task=task,
        ldscore_dir=ldscore_dir,
        annot_path=annot_path,
        sumstats_path=sumstats_path,
    )


def _run(fixture: _Fixture, tmp_path: Path) -> Path:
    scratch = tmp_path / "scratch"
    scratch.mkdir(parents=True)
    result = fixture.task.execute(
        scratch_dir=scratch, fetch=fixture.fetch, wf=make_wf()
    )
    assert isinstance(result, DirectoryAsset)
    return result.path


def test_recovers_enriched_annotation_ranking(tmp_path: Path):
    fixture = _build_fixture(tmp_path / "in", perturb_chrom=None)
    out = _run(fixture, tmp_path)
    snpvar = pl.read_parquet(out / SNPVAR_PARQUET_FILENAME)
    ann = pl.read_parquet(fixture.annot_path).rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
    joined = snpvar.join(ann, on=_KEY, how="inner")
    rho = spearmanr(joined[SNPVAR_COL], joined[ENRICHED]).statistic
    assert rho > 0.9
    # every annotation variant is scored, including those with no sumstats
    assert joined.height == snpvar.height == ann.height


def _snpvar_by_parity(out: Path, parity: int) -> pl.Series:
    return (
        pl.read_parquet(out / SNPVAR_PARQUET_FILENAME)
        .filter(pl.col(GWASLAB_CHROM_COL) % 2 == parity)
        .sort(_KEY)[SNPVAR_COL]
    )


def test_no_leakage_even_perturbation_does_not_move_even_scores(tmp_path: Path):
    # Even chromosomes are scored by tau_odd, fit on odd chromosomes only.
    # Perturbing an even chromosome's Z may move tau_even (and so the odd
    # chromosomes' snpvar) but must leave tau_odd and every even snpvar untouched.
    base = _run(_build_fixture(tmp_path / "a_in", perturb_chrom=None), tmp_path / "a")
    perturbed = _run(_build_fixture(tmp_path / "b_in", perturb_chrom=2), tmp_path / "b")
    base_diag = json.loads((base / DIAGNOSTICS_JSON_FILENAME).read_text())
    perturbed_diag = json.loads((perturbed / DIAGNOSTICS_JSON_FILENAME).read_text())
    assert base_diag["tau_odd"] == perturbed_diag["tau_odd"]
    assert base_diag["tau_even"] != perturbed_diag["tau_even"]
    assert _snpvar_by_parity(base, 0).equals(_snpvar_by_parity(perturbed, 0))
    assert not _snpvar_by_parity(base, 1).equals(_snpvar_by_parity(perturbed, 1))


def test_beta_and_se_stand_in_for_missing_z(tmp_path: Path):
    from_z = _run(_build_fixture(tmp_path / "z_in", perturb_chrom=None), tmp_path / "z")
    from_beta_se = _run(
        _build_fixture(tmp_path / "b_in", perturb_chrom=None, z_as_beta_se=True),
        tmp_path / "b",
    )
    expected = pl.read_parquet(from_z / SNPVAR_PARQUET_FILENAME)[SNPVAR_COL]
    got = pl.read_parquet(from_beta_se / SNPVAR_PARQUET_FILENAME)[SNPVAR_COL]
    assert np.allclose(got.to_numpy(), expected.to_numpy(), rtol=1e-9, atol=0)


def test_tau_weights_tables_reproduce_snpvar_exactly(tmp_path: Path):
    # Each weights table holds one half's tau in the ridge-weights schema. The
    # chromosomes that half scores must have snpvar == sum_c a_ic gamma_raw_c
    # exactly, so a contrast built from these weights explains the prior exactly.
    fixture = _build_fixture(tmp_path / "in", perturb_chrom=None)
    out = _run(fixture, tmp_path)
    snpvar = pl.read_parquet(out / SNPVAR_PARQUET_FILENAME)
    ann = pl.read_parquet(fixture.annot_path).rename(POLYFUN_TO_GWASLAB_KEY_RENAME)
    joined = snpvar.join(ann, on=_KEY, how="inner")
    # tau fit on odd chromosomes scores even ones, and vice versa.
    for filename, scored_parity in (
        (TAU_ODD_WEIGHTS_FILENAME, 0),
        (TAU_EVEN_WEIGHTS_FILENAME, 1),
    ):
        weights = pl.read_parquet(out / filename)
        assert set(weights[ANNOTATION_COL]) == set(ANNOTS + MAFBINS)
        assert weights[FAMILY_COL].null_count() == 0
        scored = joined.filter(pl.col(GWASLAB_CHROM_COL) % 2 == scored_parity)
        gamma = weights[GAMMA_RAW_COL].to_numpy()
        reconstructed = (
            scored.select(weights[ANNOTATION_COL].to_list()).to_numpy() @ gamma
        )
        assert np.allclose(reconstructed, scored[SNPVAR_COL].to_numpy(), rtol=1e-12)
