"""
System test for the SBayesRC task.

Slow-tier test: it downloads the SBayesRC example summary statistics plus the
HapMap3 EUR LD reference (~GB, via the build-system reference asset) and runs
SBayesRC through the Docker image.  Intended for scheduled / on-demand runs.
"""

import tarfile
import urllib.request
from pathlib import Path

import polars as pl

from mecfs_bio.assets.reference_data.sbayesrc.extracted.ukb_eur_hm3_ld_extracted import (
    SBAYESRC_UKB_EUR_HM3_LD_EXTRACTED,
)
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    PreformattedSBayesRCDataSource,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_task import (
    SBAYESRC_WEIGHTS_FILENAME,
    SBayesRCTask,
)
from test_mecfs_bio.system.util import StampedExternalFileTask, log_on_error

_EXAMPLE_URL = (
    "https://gctbhub.cloud.edu.au/data/SBayesRC/resources/v2.0/example/example.tar.xz"
)


def _download_example_ma(tmp_path: Path) -> Path:
    """Download and extract the SBayesRC example, returning the .ma summary file."""
    archive = tmp_path / "example.tar.xz"
    urllib.request.urlretrieve(_EXAMPLE_URL, archive)
    extract_dir = tmp_path / "example"
    extract_dir.mkdir()
    with tarfile.open(archive, "r:xz") as tar:
        tar.extractall(extract_dir)
    ma_files = sorted(extract_dir.rglob("*.ma"))
    assert ma_files, f"No .ma file found in {extract_dir}"
    return ma_files[0]


def test_sbayesrc_example(tmp_path: Path):
    """Run SBayesRC (impute -> sbayes RC) on the example data and check weights."""
    ma_path = _download_example_ma(tmp_path)

    gwas_source = PreformattedSBayesRCDataSource(
        task=StampedExternalFileTask(
            meta=GWASSummaryDataFileMeta(
                id="sbayesrc_example_ma",
                trait="sbayesrc_example",
                project="sbayesrc_example",
                sub_dir="processed",
                project_path=None,
                extension=".ma",
            ),
            external_path=ma_path,
        ),
        filename=None,
        alias="sbayesrc_example",
    )
    sbayesrc_task = SBayesRCTask.create(
        asset_id="sbayesrc_example_run",
        gwas_source=gwas_source,
        ld_reference_directory_task=SBAYESRC_UKB_EUR_HM3_LD_EXTRACTED,
        threads=4,
    )

    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"
    asset_root.mkdir(parents=True, exist_ok=True)
    with log_on_error(info_store):
        runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = runner.run([sbayesrc_task], incremental_save=True)
        assert result is not None

        out_dir = result[sbayesrc_task.asset_id].path
        weights_path = out_dir / SBAYESRC_WEIGHTS_FILENAME
        assert weights_path.is_file(), f"weights not found at {weights_path}"
        weights = pl.read_csv(weights_path, separator="\t", infer_schema_length=1000)
        # SBayesRC weights: SNP, A1, BETA, plus PIP / BETAlast columns.
        assert {"SNP", "A1", "BETA"}.issubset(set(weights.columns))
        assert weights.height > 0
