from pathlib import Path
from typing import Sequence

import pandas as pd
import pytest

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.harmonization_info import HarmonizationInfo
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.fixed_effect_meta_analysis_task import (
    META_ANALYSIS_NUM_STUDIES_COL,
    CaseControlSampleInfo,
    FixedEffectsMetaAnalysisTask,
    GwasSource,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)

# As an experiment, I asked Claude to implement unit tests
# All but one of the unit tests in this file were written by Claude
#
# Every test here drives the task through its public interface -- construct a
# FixedEffectsMetaAnalysisTask, call execute, read the output parquet -- so that the
# internal helpers stay free to change.

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

_PARQUET_READ_SPEC = DataFrameReadSpec(DataFrameParquetFormat())

_DEFAULT_SAMPLE_INFO = CaseControlSampleInfo(cases=10, controls=10)

_TRAIT = "trait"


def _harmonization_info(
    build: GenomeBuild = "19",
    ref_allele_col: str = GWASLAB_NON_EFFECT_ALLELE_COL,
) -> HarmonizationInfo:
    return HarmonizationInfo(
        build=build, ref_allele_col=ref_allele_col, pos_col=GWASLAB_POS_COL
    )


def _source(
    asset_id: str,
    sample_info: CaseControlSampleInfo = _DEFAULT_SAMPLE_INFO,
    harmonization_info: HarmonizationInfo | None = _harmonization_info(),
) -> GwasSource:
    """A source whose metadata records harmonization to a reference genome."""
    return GwasSource(
        task=FakeTask(
            FilteredGWASDataMeta(
                id=AssetId(asset_id),
                trait=_TRAIT,
                project="project",
                sub_dir="processed",
                read_spec=_PARQUET_READ_SPEC,
                harmonization_info=harmonization_info,
            )
        ),
        sample_info=sample_info,
    )


def _run_task(
    tmp_path: Path,
    dfs: Sequence[pd.DataFrame],
    sample_infos: Sequence[CaseControlSampleInfo] | None = None,
) -> pd.DataFrame:
    """Write each study to parquet, run a meta-analysis over them, return the output."""
    if sample_infos is None:
        sample_infos = [_DEFAULT_SAMPLE_INFO] * len(dfs)
    assert len(sample_infos) == len(dfs)

    paths: dict[str, Path] = {}
    sources = []
    for i, (df, sample_info) in enumerate(zip(dfs, sample_infos)):
        asset_id = f"s{i}"
        path = tmp_path / f"{asset_id}.parquet"
        df.to_parquet(path)
        paths[asset_id] = path
        sources.append(_source(asset_id, sample_info=sample_info))

    scratch = tmp_path / "scratch"
    scratch.mkdir(exist_ok=True)
    task = FixedEffectsMetaAnalysisTask(meta=SimpleFileMeta("meta"), sources=sources)

    def fetch(asset_id: AssetId) -> Asset:
        return FileAsset(paths[str(asset_id)])

    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
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


def _row_at(df: pd.DataFrame, pos: int) -> pd.Series:
    """The single output row for a given position."""
    matching = df[df[GWASLAB_POS_COL] == pos]
    assert len(matching) == 1
    return matching.iloc[0]


# ---------------------------------------------------------------------------
# Overall behaviour
# ---------------------------------------------------------------------------


def test_fixed_effect_meta_analysis_task(tmp_path: Path) -> None:
    """
    Run a meta-analysis in a simple synthetic case
    """
    df_1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=1, se=0.5)])
    df_2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=1, se=0.25)])
    df_3 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=10, se=10)])
    result = _run_task(tmp_path, [df_1, df_2, df_3])
    assert result[GWASLAB_SE_COL].item() == pytest.approx(0.2235, abs=0.01)


def test_requires_at_least_two_sources() -> None:
    """Constructing a task with a single source raises AssertionError."""
    with pytest.raises(AssertionError):
        FixedEffectsMetaAnalysisTask(
            meta=SimpleFileMeta("meta"),
            sources=[_source("s0")],
        )


# ---------------------------------------------------------------------------
# Beta weighting
# ---------------------------------------------------------------------------


def test_beta_equal_se_is_average(tmp_path: Path) -> None:
    """When both studies have equal SE, the meta-beta is their simple average."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.2, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.8, se=0.1)])
    result = _run_task(tmp_path, [df1, df2])
    assert result[GWASLAB_BETA_COL].item() == pytest.approx(0.5, abs=1e-6)
    # SE_meta = SE / sqrt(2) when SE1 == SE2
    assert result[GWASLAB_SE_COL].item() == pytest.approx(0.1 / (2**0.5), abs=1e-6)


def test_beta_unequal_se_lower_se_dominates(tmp_path: Path) -> None:
    """The study with lower SE should pull the meta-beta toward its own beta."""
    # Study 1: beta=0.0, SE=0.1  (precise, should dominate)
    # Study 2: beta=1.0, SE=0.5  (imprecise)
    # Numerator = 0/0.01 + 1/0.25 = 4
    # Denominator = 1/0.01 + 1/0.25 = 104
    # Beta = 4/104 ≈ 0.03846
    # SE   = sqrt(1/104) ≈ 0.09806
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.0, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=1.0, se=0.5)])
    result = _run_task(tmp_path, [df1, df2])
    assert result[GWASLAB_BETA_COL].item() == pytest.approx(4 / 104, abs=1e-4)
    assert result[GWASLAB_SE_COL].item() == pytest.approx((1 / 104) ** 0.5, abs=1e-4)


# ---------------------------------------------------------------------------
# Allele coding: sources are reference-harmonized, so alleles are never flipped
# ---------------------------------------------------------------------------


def test_output_keeps_input_allele_coding(tmp_path: Path) -> None:
    """Output EA/NEA and the sign of BETA are those of the harmonized inputs."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "C", "A", beta=0.5, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "C", "A", beta=0.3, se=0.1)])
    row = _row_at(_run_task(tmp_path, [df1, df2]), 10)
    assert row[GWASLAB_EFFECT_ALLELE_COL] == "C"
    assert row[GWASLAB_NON_EFFECT_ALLELE_COL] == "A"
    assert row[GWASLAB_BETA_COL] == pytest.approx(0.4, abs=1e-6)


def test_mirrored_indels_are_not_combined(tmp_path: Path) -> None:
    """An insertion and the mirrored deletion at one position are distinct variants.

    With reference GTCA at position 10, T->TCA (ref T) and TCA->T (ref TCA) are both
    consistent with the reference, and they share an unordered allele pair. Matching on
    that pair would combine the insertion's effect with the negated deletion's effect.
    """
    df1 = pd.DataFrame([_minimal_row(1, 10, "TCA", "T", beta=0.5, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "T", "TCA", beta=0.5, se=0.1)])
    result = _run_task(tmp_path, [df1, df2])
    assert len(result) == 2
    assert set(result[META_ANALYSIS_NUM_STUDIES_COL]) == {1}
    assert list(result[GWASLAB_BETA_COL]) == pytest.approx([0.5, 0.5], abs=1e-6)


# ---------------------------------------------------------------------------
# Variant set: union semantics
# ---------------------------------------------------------------------------


def test_partial_overlap_returns_union(tmp_path: Path) -> None:
    """Variants present in either study appear in the output."""
    df1 = pd.DataFrame(
        [
            _minimal_row(1, 10, "A", "C", beta=0.3, se=0.1),
            _minimal_row(1, 20, "G", "T", beta=0.5, se=0.2),
        ]
    )
    # Study 2 only has the variant at pos 10
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_task(tmp_path, [df1, df2])
    assert len(result) == 2
    # The study-1-only variant passes through untouched.
    row = _row_at(result, 20)
    assert row[GWASLAB_BETA_COL] == pytest.approx(0.5, abs=1e-6)
    assert row[GWASLAB_SE_COL] == pytest.approx(0.2, abs=1e-6)


def test_disjoint_studies_produce_all_variants(tmp_path: Path) -> None:
    """Studies with no shared variants contribute all of their variants to the union."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 20, "G", "T", beta=0.5, se=0.2)])
    result = _run_task(tmp_path, [df1, df2])
    assert len(result) == 2
    assert set(result[META_ANALYSIS_NUM_STUDIES_COL]) == {1}


def test_three_sources_with_mixed_membership(tmp_path: Path) -> None:
    """Union across three studies, with each variant backed by a different subset."""
    # pos 10: all three; pos 20: sources 0 and 2; pos 30: source 2 only.
    df1 = pd.DataFrame(
        [
            _minimal_row(1, 10, "A", "C", beta=0.2, se=0.1),
            _minimal_row(1, 20, "G", "T", beta=0.4, se=0.1),
        ]
    )
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.2, se=0.1)])
    df3 = pd.DataFrame(
        [
            _minimal_row(1, 10, "A", "C", beta=0.2, se=0.1),
            _minimal_row(1, 20, "G", "T", beta=0.4, se=0.1),
            _minimal_row(1, 30, "C", "A", beta=0.9, se=0.3),
        ]
    )
    result = _run_task(tmp_path, [df1, df2, df3])
    assert len(result) == 3

    all_three = _row_at(result, 10)
    assert all_three[META_ANALYSIS_NUM_STUDIES_COL] == 3
    assert all_three[GWASLAB_BETA_COL] == pytest.approx(0.2, abs=1e-6)
    assert all_three[GWASLAB_SE_COL] == pytest.approx(0.1 / (3**0.5), abs=1e-6)

    first_and_last = _row_at(result, 20)
    assert first_and_last[META_ANALYSIS_NUM_STUDIES_COL] == 2
    assert first_and_last[GWASLAB_EFFECT_ALLELE_COL] == "G"
    assert first_and_last[GWASLAB_NON_EFFECT_ALLELE_COL] == "T"
    assert first_and_last[GWASLAB_BETA_COL] == pytest.approx(0.4, abs=1e-6)

    last_only = _row_at(result, 30)
    assert last_only[META_ANALYSIS_NUM_STUDIES_COL] == 1
    assert last_only[GWASLAB_EFFECT_ALLELE_COL] == "C"
    assert last_only[GWASLAB_BETA_COL] == pytest.approx(0.9, abs=1e-6)


def test_output_is_sorted_by_chromosome_and_position(tmp_path: Path) -> None:
    """Output row order is deterministic, so the build cache is stable across rebuilds."""
    df1 = pd.DataFrame(
        [
            _minimal_row(2, 5, "A", "C", beta=0.3, se=0.1),
            _minimal_row(1, 30, "A", "C", beta=0.3, se=0.1),
        ]
    )
    df2 = pd.DataFrame(
        [
            _minimal_row(1, 10, "G", "T", beta=0.4, se=0.15),
            _minimal_row(2, 1, "G", "T", beta=0.4, se=0.15),
        ]
    )
    result = _run_task(tmp_path, [df1, df2])
    keys = list(zip(result[GWASLAB_CHROM_COL], result[GWASLAB_POS_COL]))
    assert keys == sorted(keys)


# ---------------------------------------------------------------------------
# Per-variant effective sample size
# ---------------------------------------------------------------------------


def test_sample_size_sums_over_all_contributing_studies(tmp_path: Path) -> None:
    """N equals the sum of effective sample sizes of the studies backing the variant."""
    # cases=100, controls=400  →  ESS = 320
    # cases=50,  controls=50   →  ESS = 100
    # total = 420
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_task(
        tmp_path,
        [df1, df2],
        [
            CaseControlSampleInfo(cases=100, controls=400),
            CaseControlSampleInfo(cases=50, controls=50),
        ],
    )
    assert result[GWASLAB_SAMPLE_SIZE_COLUMN].item() == 420
    assert result[META_ANALYSIS_NUM_STUDIES_COL].item() == 2


def test_sample_size_counts_only_contributing_studies(tmp_path: Path) -> None:
    """A variant missing from a study does not inherit that study's sample size."""
    df1 = pd.DataFrame(
        [
            _minimal_row(1, 10, "A", "C", beta=0.3, se=0.1),
            _minimal_row(1, 20, "G", "T", beta=0.3, se=0.1),
        ]
    )
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_task(
        tmp_path,
        [df1, df2],
        [
            CaseControlSampleInfo(cases=100, controls=400),  # ESS 320
            CaseControlSampleInfo(cases=50, controls=50),  # ESS 100
        ],
    )
    shared = _row_at(result, 10)
    assert shared[GWASLAB_SAMPLE_SIZE_COLUMN] == 420
    assert shared[META_ANALYSIS_NUM_STUDIES_COL] == 2

    study_one_only = _row_at(result, 20)
    assert study_one_only[GWASLAB_SAMPLE_SIZE_COLUMN] == 320
    assert study_one_only[META_ANALYSIS_NUM_STUDIES_COL] == 1


@pytest.mark.parametrize(
    ("cases", "controls", "expected_n"),
    [
        (100, 100, 200),  # equal split: ESS = total N
        (100, 400, 320),  # 4 / (1/100 + 1/400) = 320
        (1000, 9000, 3600),  # 4 / (1/1000 + 1/9000) = 3600
    ],
)
def test_effective_sample_size_of_a_single_study(
    tmp_path: Path, cases: int, controls: int, expected_n: int
) -> None:
    """The effective sample size of a case/control study, seen through a variant only
    that study has."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 20, "G", "T", beta=0.4, se=0.15)])
    result = _run_task(
        tmp_path,
        [df1, df2],
        [
            CaseControlSampleInfo(cases=cases, controls=controls),
            CaseControlSampleInfo(cases=7, controls=7),
        ],
    )
    assert _row_at(result, 10)[GWASLAB_SAMPLE_SIZE_COLUMN] == expected_n


# ---------------------------------------------------------------------------
# rsID
# ---------------------------------------------------------------------------


def _row_with_rsid(
    chrom: int, pos: int, ea: str, nea: str, beta: float, se: float, rsid: str
) -> dict:
    return _minimal_row(chrom, pos, ea, nea, beta, se) | {GWASLAB_RSID_COL: rsid}


def test_rsid_from_study1_is_preserved(tmp_path: Path) -> None:
    """rsID present in study 1 flows through to the output."""
    df1 = pd.DataFrame(
        [_row_with_rsid(1, 10, "A", "C", beta=0.3, se=0.1, rsid="rs12345")]
    )
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_task(tmp_path, [df1, df2])
    assert GWASLAB_RSID_COL in result.columns
    assert result[GWASLAB_RSID_COL].item() == "rs12345"


def test_rsid_recovered_from_later_study(tmp_path: Path) -> None:
    """A variant absent from study 1 still gets its rsID from a study that has it."""
    df1 = pd.DataFrame(
        [_row_with_rsid(1, 10, "A", "C", beta=0.3, se=0.1, rsid="rs12345")]
    )
    df2 = pd.DataFrame(
        [_row_with_rsid(1, 20, "G", "T", beta=0.4, se=0.15, rsid="rs67890")]
    )
    result = _run_task(tmp_path, [df1, df2])
    assert _row_at(result, 20)[GWASLAB_RSID_COL] == "rs67890"


def test_no_rsid_column_when_no_study_has_it(tmp_path: Path) -> None:
    """When no study carries an rsID column the output has none either."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    result = _run_task(tmp_path, [df1, df2])
    assert GWASLAB_RSID_COL not in result.columns


# ---------------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------------


def test_duplicate_variant_within_a_study_is_rejected(tmp_path: Path) -> None:
    """The same variant twice in one study would double-count its weight."""
    df1 = pd.DataFrame(
        [
            _minimal_row(1, 10, "A", "C", beta=0.3, se=0.1),
            _minimal_row(1, 10, "A", "C", beta=0.5, se=0.1),
        ]
    )
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    with pytest.raises(AssertionError):
        _run_task(tmp_path, [df1, df2])


def test_zero_standard_error_is_rejected(tmp_path: Path) -> None:
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.0)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    with pytest.raises(AssertionError):
        _run_task(tmp_path, [df1, df2])


def test_null_standard_error_is_rejected(tmp_path: Path) -> None:
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df1[GWASLAB_SE_COL] = pd.array([None], dtype="Float64")
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    with pytest.raises(AssertionError):
        _run_task(tmp_path, [df1, df2])


def test_null_beta_is_rejected(tmp_path: Path) -> None:
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.3, se=0.1)])
    df1[GWASLAB_BETA_COL] = pd.array([None], dtype="Float64")
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    with pytest.raises(AssertionError):
        _run_task(tmp_path, [df1, df2])


def test_effect_allele_equal_to_non_effect_allele_is_rejected(tmp_path: Path) -> None:
    """A variant whose two alleles are identical has no well-defined orientation."""
    df1 = pd.DataFrame([_minimal_row(1, 10, "A", "A", beta=0.3, se=0.1)])
    df2 = pd.DataFrame([_minimal_row(1, 10, "A", "C", beta=0.4, se=0.15)])
    with pytest.raises(AssertionError):
        _run_task(tmp_path, [df1, df2])


# ---------------------------------------------------------------------------
# Harmonization requirements
# ---------------------------------------------------------------------------


def test_unharmonized_source_is_rejected() -> None:
    with pytest.raises(AssertionError):
        FixedEffectsMetaAnalysisTask.create(
            asset_id="meta",
            meta_analysis_name="meta",
            sources=[_source("s0"), _source("s1", harmonization_info=None)],
        )


def test_sources_harmonized_to_different_builds_are_rejected() -> None:
    with pytest.raises(AssertionError):
        FixedEffectsMetaAnalysisTask.create(
            asset_id="meta",
            meta_analysis_name="meta",
            sources=[
                _source("s0", harmonization_info=_harmonization_info(build="19")),
                _source("s1", harmonization_info=_harmonization_info(build="38")),
            ],
        )


def test_source_with_reference_allele_in_effect_allele_column_is_rejected() -> None:
    with pytest.raises(AssertionError):
        FixedEffectsMetaAnalysisTask.create(
            asset_id="meta",
            meta_analysis_name="meta",
            sources=[
                _source("s0"),
                _source(
                    "s1",
                    harmonization_info=_harmonization_info(
                        ref_allele_col=GWASLAB_EFFECT_ALLELE_COL
                    ),
                ),
            ],
        )


def test_output_is_harmonized_to_the_sources_build() -> None:
    """The output can feed tasks that require harmonized input, such as fine-mapping."""
    task = FixedEffectsMetaAnalysisTask.create(
        asset_id="meta",
        meta_analysis_name="meta",
        sources=[
            _source("s0", harmonization_info=_harmonization_info(build="38")),
            _source("s1", harmonization_info=_harmonization_info(build="38")),
        ],
    )
    assert isinstance(task.meta, FilteredGWASDataMeta)
    assert task.meta.harmonization_info == _harmonization_info(build="38")
