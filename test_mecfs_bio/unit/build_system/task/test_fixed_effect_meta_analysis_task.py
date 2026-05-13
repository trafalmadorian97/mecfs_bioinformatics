from pathlib import Path

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fixed_effect_meta_analysis_task import (
    CaseControlSampleInfo,
    FixedEffectsMetaAnalysisTask,
    GwasSource,
)
from mecfs_bio.build_system.wf.base_wf import SimpleWF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECTIVE_SAMPLE_SIZE,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL,
)

# As an experiment, I asked Claude to implement unit tests
# All but one of the unit tests in this file were written by Claude

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PARQUET_READ_SPEC = DataFrameReadSpec(DataFrameParquetFormat())


def _make_source(asset_id: str, cases: int, controls: int) -> GwasSource:
    return GwasSource(
        task=FakeTask(SimpleFileMeta(asset_id, read_spec=_PARQUET_READ_SPEC)),
        sample_info=CaseControlSampleInfo(cases=cases, controls=controls),
    )


def _run_two_source_task(
    tmp_path: Path,
    df1: pd.DataFrame,
    df2: pd.DataFrame,
    cases1: int = 10,
    controls1: int = 10,
    cases2: int = 10,
    controls2: int = 10,
) -> pd.DataFrame:
    """Write two parquet files and run a two-source meta-analysis, returning the result."""
    p1 = tmp_path / "s1.parquet"
    p2 = tmp_path / "s2.parquet"
    df1.to_parquet(p1)
    df2.to_parquet(p2)
    scratch = tmp_path / "scratch"
    scratch.mkdir(exist_ok=True)

    task = FixedEffectsMetaAnalysisTask(
        meta=SimpleFileMeta("meta"),
        sources=[
            _make_source("s1", cases=cases1, controls=controls1),
            _make_source("s2", cases=cases2, controls=controls2),
        ],
    )

    paths = {"s1": p1, "s2": p2}

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(paths[str(asset_id)])

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    return pd.read_parquet(result.path)


def _minimal_row(
    chrom: int, pos: int, ea: str, nea: str, beta: float, se: float
) -> dict:
    return {
        GWASLAB_CHROM_COL: chrom,
        GWASLAB_POS_COL: pos,
        GWASLAB_EFFECT_ALLELE_COL: ea,
        GWASLAB_NON_EFFECT_ALLELE_COL: nea,
        GWASLAB_BETA_COL: beta,
        GWASLAB_SE_COL: se,
    }


def test_fixed_effect_meta_analysis_task(tmp_path: Path):
    """
    Run a meta-analysis in a simple synthetic case
    """
    df_1 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [1],
            GWASLAB_SE_COL: [0.5],
        }
    )

    df_2 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [1],
            GWASLAB_SE_COL: [0.25],
        }
    )
    df_3 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [10],
            GWASLAB_SE_COL: [10],
        }
    )
    df_1_path = tmp_path / "df_1.parquet"
    df_1.to_parquet(df_1_path)
    df_2_path = tmp_path / "df_2.parquet"
    df_2.to_parquet(df_2_path)
    df_3_path = tmp_path / "df_3.parquet"
    df_3.to_parquet(df_3_path)
    scratch = tmp_path / "scratch"
    scratch.mkdir()
    task = FixedEffectsMetaAnalysisTask(
        meta=SimpleFileMeta("meta_analysis"),
        sources=[
            GwasSource(
                task=FakeTask(
                    SimpleFileMeta(
                        "df1", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                sample_info=CaseControlSampleInfo(controls=10, cases=10),
            ),
            GwasSource(
                task=FakeTask(
                    SimpleFileMeta(
                        "df2", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                sample_info=CaseControlSampleInfo(controls=10, cases=10),
            ),
            GwasSource(
                task=FakeTask(
                    SimpleFileMeta(
                        "df3", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
                    )
                ),
                sample_info=CaseControlSampleInfo(controls=10, cases=10),
            ),
        ],
    )

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "df1":
            return FileAsset(df_1_path)
        if asset_id == "df2":
            return FileAsset(df_2_path)
        if asset_id == "df3":
            return FileAsset(df_3_path)
        raise ValueError("unknown asset id")

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=SimpleWF())
    assert isinstance(result, FileAsset)
    df = pd.read_parquet(result.path)
    assert df["SE"].item() == pytest.approx(0.2235, abs=0.01)


# ---------------------------------------------------------------------------
# Effective sample size math
# ---------------------------------------------------------------------------


def test_case_control_effective_sample_size_equal_split():
    """Equal cases/controls: ESS = N (equal to total N)."""
    info = CaseControlSampleInfo(cases=100, controls=100)
    assert info.effective_sample_size() == 200


def test_case_control_effective_sample_size_unequal_split():
    """ESS < N when cases != controls (harmonic mean penalty)."""
    # 4 / (1/100 + 1/400) = 4 / 0.0125 = 320
    info = CaseControlSampleInfo(cases=100, controls=400)
    assert info.effective_sample_size() == 320


def test_case_control_effective_sample_size_extreme_imbalance():
    # 4 / (1/1000 + 1/9000) = 4 / (10/9000) = 3600
    info = CaseControlSampleInfo(cases=1000, controls=9000)
    assert info.effective_sample_size() == 3600


# ---------------------------------------------------------------------------
# Require at least two sources
# ---------------------------------------------------------------------------


def test_requires_at_least_two_sources():
    """Constructing a task with a single source raises AssertionError."""
    with pytest.raises(AssertionError):
        FixedEffectsMetaAnalysisTask(
            meta=SimpleFileMeta("meta"),
            sources=[_make_source("s1", cases=10, controls=10)],
        )


# ---------------------------------------------------------------------------
# Beta weighting
# ---------------------------------------------------------------------------


def test_beta_equal_se_is_average(tmp_path: Path):
    """When both studies have equal SE, the meta-beta is their simple average."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.2, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.8, se=0.1)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert result[GWASLAB_BETA_COL].item() == pytest.approx(0.5, abs=1e-6)
    # SE_meta = SE / sqrt(2) when SE1 == SE2
    assert result[GWASLAB_SE_COL].item() == pytest.approx(0.1 / (2**0.5), abs=1e-6)


def test_beta_unequal_se_lower_se_dominates(tmp_path: Path):
    """The study with lower SE should pull the meta-beta toward its own beta."""
    # Study 1: beta=0.0, SE=0.1  (precise, should dominate)
    # Study 2: beta=1.0, SE=0.5  (imprecise)
    # Numerator = 0/0.01 + 1/0.25 = 4
    # Denominator = 1/0.01 + 1/0.25 = 104
    # Beta = 4/104 ≈ 0.03846
    # SE   = sqrt(1/104) ≈ 0.09806
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.0, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=1.0, se=0.5)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert result[GWASLAB_BETA_COL].item() == pytest.approx(4 / 104, abs=1e-4)
    assert result[GWASLAB_SE_COL].item() == pytest.approx((1 / 104) ** 0.5, abs=1e-4)


# ---------------------------------------------------------------------------
# Allele flip harmonisation
# ---------------------------------------------------------------------------


def test_allele_flip_harmonises_beta(tmp_path: Path):
    """Study 2 reported with swapped alleles (opposite strand coding) is correctly harmonised."""
    # Study 1: EA=A, NEA=C, BETA=+0.5  → A increases trait
    # Study 2: EA=C, NEA=A, BETA=-0.5  → same biology, opposite coding
    # After flip study 2 contributes BETA=+0.5 with SE=0.1
    # Meta-beta should equal 0.5, meta-SE = 0.1/sqrt(2)
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.5, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "C", "A", beta=-0.5, se=0.1)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert len(result) == 1
    assert result[GWASLAB_BETA_COL].item() == pytest.approx(0.5, abs=1e-6)
    assert result[GWASLAB_SE_COL].item() == pytest.approx(0.1 / (2**0.5), abs=1e-6)


def test_allele_flip_sets_flipped_flag(tmp_path: Path):
    """The flipped_1 column is True for a variant that needed harmonisation."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.5, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "C", "A", beta=-0.5, se=0.1)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert result["flipped_1"].item() is True


def test_no_flip_sets_flipped_flag_false(tmp_path: Path):
    """The flipped_1 column is False when alleles already match without flipping."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.5, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.5, se=0.1)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert result["flipped_1"].item() is False


# ---------------------------------------------------------------------------
# Variant set: intersection semantics
# ---------------------------------------------------------------------------


def test_partial_overlap_returns_intersection(tmp_path: Path):
    """Only variants shared by both studies appear in the output."""
    df1 = pd.DataFrame(
        [
            _minimal_row(1, 10, "A", "C", beta=0.3, se=0.1),
            _minimal_row(1, 20, "G", "T", beta=0.5, se=0.2),
        ]
    )
    # Study 2 only has the variant at pos 10
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert len(result) == 1
    assert result[GWASLAB_POS_COL].item() == 10


def test_empty_intersection_produces_empty_output(tmp_path: Path):
    """Studies with no shared variants produce a zero-row result (no error)."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 20, "G", "T", beta=0.5, se=0.2)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert len(result) == 0


# ---------------------------------------------------------------------------
# Additional columns
# ---------------------------------------------------------------------------


def test_n_eff_column_is_correct(tmp_path: Path):
    """N_EFF equals the sum of effective sample sizes across all sources."""
    # cases=100, controls=400  →  ESS = 320
    # cases=50,  controls=50   →  ESS = 100
    # total = 420
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_two_source_task(
        tmp_path, df1, df2, cases1=100, controls1=400, cases2=50, controls2=50
    )
    assert result[GWASLAB_EFFECTIVE_SAMPLE_SIZE].item() == 420


def test_rsid_from_study1_is_preserved(tmp_path: Path):
    """rsID present in study 1 flows through to the output."""
    df1 = pd.DataFrame(
        {
            GWASLAB_CHROM_COL: [1],
            GWASLAB_POS_COL: [10],
            GWASLAB_EFFECT_ALLELE_COL: ["A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["C"],
            GWASLAB_BETA_COL: [0.3],
            GWASLAB_SE_COL: [0.1],
            GWASLAB_RSID_COL: ["rs12345"],
        }
    )
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert GWASLAB_RSID_COL in result.columns
    assert result[GWASLAB_RSID_COL].item() == "rs12345"


def test_no_rsid_column_when_study1_lacks_it(tmp_path: Path):
    """When study 1 has no rsID column the output also has no rsID column."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_two_source_task(tmp_path, df1, df2)
    assert GWASLAB_RSID_COL not in result.columns
