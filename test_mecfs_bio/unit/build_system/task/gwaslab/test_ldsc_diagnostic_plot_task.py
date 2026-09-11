import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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
    QuantPhenotype,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_snp_heritability_by_ldsc_task import (
    SNPHeritabilityByLDSCTask,
)
from mecfs_bio.build_system.task.gwaslab.ldsc_diagnostic_plot_task import (
    LdscDiagnosticPlotConfig,
    LdscDiagnosticPlotTask,
    LdscFit,
    compute_chi2,
    estimate_observed_fit,
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


def test_estimate_observed_fit_recovers_a_noise_free_line():
    # Points lying exactly on E[chi2] = 1 + (N*h2/M)*ld must be recovered by the fit,
    # regardless of the weighting: intercept 1.0 and, with slope 0.2 = 1000*0.2/1000, h2 0.2.
    n, m = 1000.0, 1000.0
    ld = np.linspace(1.0, 300.0, 500)
    chi2 = 1.0 + (n * 0.2 / m) * ld
    fit = estimate_observed_fit(chi2=chi2, ld=ld, n=n, m=m)
    assert fit.intercept == pytest.approx(1.0, abs=1e-6)
    assert fit.h2_obs == pytest.approx(0.2, abs=1e-6)


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
        sumstats_reader=lambda _asset: sumstats_df,
        fit_estimator=lambda chi2, ld, n, m: LdscFit(intercept=1.0, h2_obs=0.2),
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
