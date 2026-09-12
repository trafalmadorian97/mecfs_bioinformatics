import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from gwaslab.extension.ldsc.ldsc_regressions import h2_obs_to_liab

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    SNPHeritabilityByLDSCTask,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic import bin_by_ld_score
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import (
    LdscDiagnosticPlotConfig,
    LdscDiagnosticPlotTask,
    LdscFit,
    build_diagnostic_figure,
    compute_chi2,
    gwaslab_observed_fit,
    liability_h2,
    merge_chi2_with_ld_scores,
)
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SAMPLE_SIZE_COLUMN,
    GWASLAB_SE_COL,
)
from mecfs_bio.constants.ldsc_constants import LDSC_Z_COL


def test_compute_chi2_prefers_z_column():
    df = pd.DataFrame({LDSC_Z_COL: [2.0, -3.0]})
    assert compute_chi2(df).tolist() == [4.0, 9.0]


def test_compute_chi2_falls_back_to_beta_over_se():
    df = pd.DataFrame({GWASLAB_BETA_COL: [2.0, -3.0], GWASLAB_SE_COL: [1.0, 1.0]})
    assert compute_chi2(df).tolist() == [4.0, 9.0]


def test_merge_matches_on_rsid_and_aligns_chi2_with_ld():
    sumstats = pd.DataFrame(
        {GWASLAB_RSID_COL: ["rs1", "rs2", "rs3"], LDSC_Z_COL: [1.0, 2.0, 3.0]}
    )
    # rs4 has no sumstats row and rs1 no LD score, so only rs2/rs3 survive the inner join.
    ld_scores = pd.DataFrame({"SNP": ["rs2", "rs3", "rs4"], "L2": [10.0, 20.0, 30.0]})
    merged = merge_chi2_with_ld_scores(sumstats, ld_scores)
    assert merged.chi2.tolist() == [4.0, 9.0]
    assert merged.ld.tolist() == [10.0, 20.0]


class _FakeSumstats:
    """One stand-in satisfying LdscSumstats: holds a variant table (for the binning path), records
    the estimate call, and serves a canned ldsc_h2 summary -- so both the fit's column extraction
    and the task's execute can be tested without running LDSC or reading an LD reference."""

    def __init__(
        self,
        data: pd.DataFrame | None = None,
        ldsc_h2_summary: pd.DataFrame | None = None,
        build: str = "19",
    ):
        self.data = data if data is not None else pd.DataFrame()
        self.meta = {"gwaslab": {"genome_build": build}}
        self._summary = ldsc_h2_summary
        self.ldsc_h2: pd.DataFrame | None = None
        self.estimate_call: dict | None = None

    def infer_build(self) -> None:
        pass

    def estimate_h2_by_ldsc(
        self,
        *,
        ref_ld_chr: str,
        w_ld_chr: str,
        samp_prev: float | None,
        pop_prev: float | None,
    ) -> None:
        self.estimate_call = dict(
            ref_ld_chr=ref_ld_chr, samp_prev=samp_prev, pop_prev=pop_prev
        )
        self.ldsc_h2 = self._summary


def test_gwaslab_observed_fit_reads_intercept_and_observed_h2():
    # gwaslab's parse_ldsc_summary yields string-valued Intercept and h2_obs columns.
    fake = _FakeSumstats(
        ldsc_h2_summary=pd.DataFrame({"Intercept": ["0.93"], "h2_obs": ["0.2"]}),
        build="19",
    )
    fit = gwaslab_observed_fit(fake, ref_ld_chr="/ref/LDscore.@", build="19")
    assert fit.intercept == pytest.approx(0.93)
    assert fit.h2_obs == pytest.approx(0.2)


def test_gwaslab_observed_fit_requests_observed_scale():
    # Observed scale means None (not NaN) prevalences: gwaslab labels the result "Liability" and
    # renames the h2 column whenever both prevalences are not None, so None is what keeps the plain
    # h2_obs column -- the scale in which the fitted slope is reconstructable.
    fake = _FakeSumstats(
        ldsc_h2_summary=pd.DataFrame({"Intercept": ["1.0"], "h2_obs": ["0.1"]})
    )
    gwaslab_observed_fit(fake, ref_ld_chr="ref", build="19")
    assert fake.estimate_call is not None
    assert fake.estimate_call["samp_prev"] is None
    assert fake.estimate_call["pop_prev"] is None


def test_liability_h2_is_none_for_a_quantitative_phenotype():
    assert liability_h2(0.1, QuantPhenotype(total_sample_size=1000)) is None


def test_liability_h2_applies_gwaslabs_conversion_with_sample_and_population_prevalence():
    # P is the sample prevalence, K the population prevalence. Choosing P != K makes the
    # conversion asymmetric, so a swapped-argument regression would not match.
    pheno = BinaryPhenotypeSampleInfo(
        sample_prevalence=0.3, estimated_population_prevalence=0.05
    )
    expected = h2_obs_to_liab(0.1, 0.3, 0.05)
    assert liability_h2(0.1, pheno) == pytest.approx(expected)


def test_figure_annotation_includes_liability_line_only_when_given():
    bins = bin_by_ld_score(
        np.array([1.0, 2.0, 3.0, 4.0]), np.array([1.0, 2.0, 3.0, 4.0]), n_bins=2
    )
    fit = LdscFit(intercept=1.0, h2_obs=0.2)
    config = LdscDiagnosticPlotConfig(n_bins=2)
    with_liab = build_diagnostic_figure(
        bins=bins, fit=fit, n=1000.0, m=500.0, config=config, h2_liability=0.45
    )
    without_liab = build_diagnostic_figure(
        bins=bins, fit=fit, n=1000.0, m=500.0, config=config, h2_liability=None
    )
    assert "liability" in with_liab.layout.annotations[0].text.lower()
    assert "liability" not in without_liab.layout.annotations[0].text.lower()


def test_figure_anchors_a_trace_to_the_secondary_axis():
    # A secondary y-axis only renders when a trace references it, so guard that one does -- without
    # it the chi^2 * M / N axis silently disappears.
    bins = bin_by_ld_score(
        np.array([1.0, 2.0, 3.0, 4.0]), np.array([1.0, 2.0, 3.0, 4.0]), n_bins=2
    )
    fig = build_diagnostic_figure(
        bins=bins,
        fit=LdscFit(intercept=1.0, h2_obs=0.2),
        n=1000.0,
        m=500.0,
        config=LdscDiagnosticPlotConfig(n_bins=2),
    )
    assert fig.layout.yaxis2.overlaying == "y"
    assert any(trace.yaxis == "y2" for trace in fig.data)


def _write_ld_dir(ld_dir: Path, rsids: list[str], l2: list[float]) -> None:
    ld_dir.mkdir(parents=True, exist_ok=True)
    lines = ["CHR\tSNP\tBP\tL2\n"]
    for i, (rsid, score) in enumerate(zip(rsids, l2)):
        lines.append(f"1\t{rsid}\t{100 + i}\t{score}\n")
    with gzip.open(ld_dir / "LDscore.1.l2.ldscore.gz", "wt") as f:
        f.writelines(lines)
    (ld_dir / "LDscore.1.l2.M_5_50").write_text("1000\n")


def test_task_execute_writes_a_plot(tmp_path: Path):
    rsids = [f"rs{i}" for i in range(9)]
    ld_dir = tmp_path / "ldref"
    _write_ld_dir(ld_dir, rsids, l2=[float(i + 1) for i in range(9)])

    sumstats_df = pd.DataFrame(
        {
            GWASLAB_RSID_COL: rsids,
            GWASLAB_BETA_COL: [0.1 * (i + 1) for i in range(9)],
            GWASLAB_SE_COL: [1.0] * 9,
            GWASLAB_SAMPLE_SIZE_COLUMN: [1000] * 9,
        }
    )

    ldsc_task = SNPHeritabilityByLDSCTask.create(
        asset_id="h2",
        source_sumstats_task=FakeTask(
            GWASLabSumStatsMeta(id=AssetId("sumstats"), trait="trait", project="proj")
        ),
        ld_ref_task=FakeTask(SimpleDirectoryMeta("ldref")),
        phenotype_info=QuantPhenotype(total_sample_size=1000),
        pipe=IdentityPipe(),
        build="19",
    )
    task = LdscDiagnosticPlotTask.create(
        asset_id="ldsc_diagnostic",
        ldsc_task=ldsc_task,
        config=LdscDiagnosticPlotConfig(n_bins=3),
        sumstats_reader=lambda _asset: _FakeSumstats(data=sumstats_df),
        fit_estimator=lambda _sumstats, _ref, _build: LdscFit(
            intercept=1.0, h2_obs=0.2
        ),
    )

    # The faked sumstats_reader ignores the asset, but FileAsset still asserts the file exists.
    dummy_sumstats = tmp_path / "unused_sumstats.pkl"
    dummy_sumstats.write_bytes(b"")

    def fetch(asset_id: AssetId) -> Asset:
        if asset_id == "ldref":
            return DirectoryAsset(ld_dir)
        return FileAsset(dummy_sumstats)

    scratch = tmp_path / "scratch"
    scratch.mkdir()
    result = task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())
    assert isinstance(result, DirectoryAsset)
    assert list(result.path.glob("**/*.html"))
