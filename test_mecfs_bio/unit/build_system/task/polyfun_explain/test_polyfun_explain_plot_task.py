import shutil
from pathlib import Path

import polars as pl
import pytest
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameParquetFormat,
    DataFrameReadSpec,
)
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.polyfun_explain.polyfun_explain_plot_task import (
    PLOT_PNG_FILENAME,
    PLOT_SVG_FILENAME,
    PolyfunExplainPlotTask,
    PriorRun,
    _wrap_callout_label,
)
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    FILTERED_GWAS_FILENAME,
)
from mecfs_bio.build_system.task.susie_stacked_plot_task import (
    GENE_INFO_CHROM_COL,
    GENE_INFO_END_COL,
    GENE_INFO_NAME_COL,
    GENE_INFO_START_COL,
    GENE_INFO_STRAND_COL,
)
from mecfs_bio.build_system.wf.base_wf import make_wf
from mecfs_bio.constants.genetic_map_constants import (
    GMAP_CM_COL,
    GMAP_POS_COL,
    GMAP_RATE_COL,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL
from test_mecfs_bio.unit.build_system.task.polyfun_explain.test_polyfun_explain_contrast_task import (
    build_synthetic_explain_inputs,
)


def test_wrap_callout_label_short_unchanged_long_wrapped():
    fonts = (10.0, 9.5)
    short = "173855298:A:T (conserved ++, coding +)"
    text, _size = _wrap_callout_label(short, fonts)
    assert text == short

    long = "47731228:A:C (coding ++, ld related continuous +, open chromatin +)"
    text, _size = _wrap_callout_label(long, fonts)
    assert text == (
        "47731228:A:C\ncoding ++\nld related continuous +\nopen chromatin +"
    )


@frozen(slots=True)
class _ReferenceTracks:
    """Gene and genetic-map inputs covering the synthetic locus (POS 10-60)."""

    gene_task: FakeTask
    genetic_map_task: FakeTask
    assets: dict[str, Asset]


def _reference_tracks(tmp_path: Path) -> _ReferenceTracks:
    gene_info = pl.DataFrame(
        {
            GENE_INFO_CHROM_COL: [1],
            GENE_INFO_START_COL: [5],
            GENE_INFO_END_COL: [65],
            GENE_INFO_STRAND_COL: ["+"],
            GENE_INFO_NAME_COL: ["GENE1"],
        }
    )
    gene_path = tmp_path / "genes.parquet"
    gene_info.write_parquet(gene_path)
    # hg19 genetic map covering the locus; the recomb track reads the rate
    # column directly.
    genetic_map = pl.DataFrame(
        {
            "CHR": [1, 1, 1],
            GMAP_POS_COL: [10, 35, 60],
            GMAP_RATE_COL: [0.5, 2.0, 1.0],
            GMAP_CM_COL: [0.0, 0.5, 1.2],
        }
    )
    gmap_path = tmp_path / "genetic_map.parquet"
    genetic_map.write_parquet(gmap_path)
    return _ReferenceTracks(
        gene_task=FakeTask(
            SimpleFileMeta(
                "genes", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            )
        ),
        genetic_map_task=FakeTask(
            SimpleFileMeta(
                "genetic_map", read_spec=DataFrameReadSpec(DataFrameParquetFormat())
            )
        ),
        assets={"genes": FileAsset(gene_path), "genetic_map": FileAsset(gmap_path)},
    )


def _execute_plot(
    tmp_path: Path,
    uni_task: Task,
    prior_runs: tuple[PriorRun, ...],
    assets: dict[str, Asset],
) -> Asset:
    tracks = _reference_tracks(tmp_path)
    plot_task = PolyfunExplainPlotTask.create(
        asset_id="plot",
        susie_uniform_task=uni_task,
        prior_runs=prior_runs,
        gene_info_task=tracks.gene_task,
        genetic_map_task=tracks.genetic_map_task,
        gene_info_pipe=IdentityPipe(),
    )
    fetch_map = {**assets, **tracks.assets}

    def fetch(asset_id: AssetId) -> Asset:
        return fetch_map[str(asset_id)]

    scratch = tmp_path / "plot_scratch"
    scratch.mkdir()
    return plot_task.execute(scratch_dir=scratch, fetch=fetch, wf=make_wf())


@pytest.mark.parametrize("n_prior_runs", [1, 2])
def test_plot_writes_png_and_svg(tmp_path: Path, n_prior_runs: int):
    inputs = build_synthetic_explain_inputs(tmp_path)
    # The fixture has one prior run; drawing it under several labels exercises
    # the one-row-per-prior layout without a second synthetic run.
    prior_runs = tuple(
        PriorRun(
            label=f"prior {i}",
            susie_task=inputs.pf_task,
            contrast_task=inputs.contrast_task,
        )
        for i in range(n_prior_runs)
    )
    result = _execute_plot(
        tmp_path, inputs.uni_task, prior_runs, assets=dict(inputs.fetch_map)
    )
    assert isinstance(result, DirectoryAsset)
    png_path = result.path / PLOT_PNG_FILENAME
    svg_path = result.path / PLOT_SVG_FILENAME
    assert png_path.is_file()
    assert svg_path.is_file()
    assert png_path.stat().st_size > 0
    assert svg_path.stat().st_size > 0


def test_prior_run_on_another_chromosome_fails(tmp_path: Path):
    inputs = build_synthetic_explain_inputs(tmp_path)
    assets = dict(inputs.fetch_map)
    # Copy the prior run's directory with its variants moved to chromosome 2.
    pf_dir = assets[str(inputs.pf_task.asset_id)]
    assert isinstance(pf_dir, DirectoryAsset)
    other_dir = tmp_path / "other_chrom_prior"
    shutil.copytree(pf_dir.path, other_dir)
    gwas_path = other_dir / FILTERED_GWAS_FILENAME
    pl.read_parquet(gwas_path).with_columns(
        pl.lit(2).cast(pl.Int64).alias(GWASLAB_CHROM_COL)
    ).write_parquet(gwas_path)
    other_task = FakeTask(
        ResultDirectoryMeta(id=AssetId("other_chrom_prior"), trait="t", project="p")
    )
    assets["other_chrom_prior"] = DirectoryAsset(other_dir)
    prior_runs = (
        PriorRun(
            label="other chromosome",
            susie_task=other_task,
            contrast_task=inputs.contrast_task,
        ),
    )
    with pytest.raises(AssertionError):
        _execute_plot(tmp_path, inputs.uni_task, prior_runs, assets=assets)
