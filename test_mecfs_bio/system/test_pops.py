"""
System test for the POPs tasks: munge features, then run POPs.

Uses the small example data bundled inside the POPs source tarball (the
PASS_Schizophrenia example from the POPs tutorial): 8 raw feature files over ~18k
genes and a matching MAGMA gene-analysis result. Exercises both
PopsMungeFeatureDirectoryTask and PopsRunTask end to end and verifies that a
non-empty .preds table is produced.
"""

from pathlib import Path

import numpy as np
import pandas as pd

from mecfs_bio.assets.reference_data.pops.source.pops_source_extracted import (
    POPS_SOURCE_EXTRACTED,
)
from mecfs_bio.assets.reference_data.pops.source.pops_source_raw import (
    POPS_SOURCE_RAW,
    POPS_SOURCE_TARBALL_TOP_DIR,
)
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask
from mecfs_bio.build_system.task.pops.pops_munge_task import (
    PopsMungeFeatureDirectoryTask,
)
from mecfs_bio.build_system.task.pops.pops_run_task import PopsRunTask
from mecfs_bio.build_system.task.pops.pops_utils import (
    CONTROL_FEATURES_RELATIVE_PATH,
    GENE_ANNOT_RELATIVE_PATH,
    POPS_OUTPUT_STEM_NAME,
    POPS_PREDS_SUFFIX,
    POPS_SCRIPT_NAME,
    invoke_pops_lowmem_script,
    invoke_pops_script,
)
from test_mecfs_bio.system.util import log_on_error

# The example MAGMA scores in the tarball use this stem (not the stem our own MAGMA
# pipeline uses).
_EXAMPLE_MAGMA_STEM = "PASS_Schizophrenia"

# Sub-directories of the POPs tarball holding the tutorial example data.
_EXAMPLE_FEATURES_RAW_TASK = ExtractTarGzipTask.create(
    asset_id="pops_example_features_raw",
    source_task=POPS_SOURCE_RAW,
    sub_folder_name_inside_tar=(
        f"{POPS_SOURCE_TARBALL_TOP_DIR}/example/data/features_raw"
    ),
)
_EXAMPLE_MAGMA_SCORES_TASK = ExtractTarGzipTask.create(
    asset_id="pops_example_magma_scores",
    source_task=POPS_SOURCE_RAW,
    sub_folder_name_inside_tar=(
        f"{POPS_SOURCE_TARBALL_TOP_DIR}/example/data/magma_scores"
    ),
)


def test_pops_munge_and_run_example(tmp_path: Path):
    """Munge the example raw features and run POPs against the example MAGMA scores,
    verifying that a non-empty .preds table with a PoPS score column is written."""
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    munge_task = PopsMungeFeatureDirectoryTask.create(
        asset_id="pops_example_munge",
        pops_source_task=POPS_SOURCE_EXTRACTED,
        raw_features_task=_EXAMPLE_FEATURES_RAW_TASK,
        max_cols=500,
    )
    pops_run_task = PopsRunTask(
        meta=SimpleDirectoryMeta(id=AssetId("pops_example_run")),
        pops_source_task=POPS_SOURCE_EXTRACTED,
        munged_features_task=munge_task,
        magma_gene_analysis_task=_EXAMPLE_MAGMA_SCORES_TASK,
        magma_stem_name=_EXAMPLE_MAGMA_STEM,
    )

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run([munge_task, pops_run_task], incremental_save=True)
        assert result is not None

        munge_output = result[munge_task.asset_id]
        assert isinstance(munge_output, DirectoryAsset)
        chunk_cols = list(munge_output.path.glob("pops_features.cols.*.txt"))
        assert len(chunk_cols) > 0, "munge produced no feature chunks"
        assert (munge_output.path / "pops_features.rows.txt").exists()

        pops_output = result[pops_run_task.asset_id]
        assert isinstance(pops_output, DirectoryAsset)
        preds_path = pops_output.path / (POPS_OUTPUT_STEM_NAME + POPS_PREDS_SUFFIX)
        assert preds_path.exists(), f"POPs .preds output not found at {preds_path}"
        lines = preds_path.read_text().splitlines()
        assert len(lines) > 1, "POPs .preds should have a header and gene rows"
        header = lines[0]
        assert "ENSGID" in header
        assert "PoPS_Score" in header


def _write_chr22_magma_subset(example_magma_dir: Path, out_dir: Path) -> Path:
    """Derive a small single-chromosome (chr22, ~394 genes) MAGMA subset from the
    tarball's PASS_Schizophrenia example so the low-memory POPs fit is fast enough
    for CI (its Gram eigendecomposition is O(n_genes cubed)). Returns the subset
    MAGMA prefix. genes.out keeps its header plus chr22 rows; genes.raw keeps its two
    comment headers plus chr22 gene lines (whose covariances reference only chr22
    genes, so the block stays self-consistent)."""
    out_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = out_dir / "chr22"
    genes_out = (
        (example_magma_dir / f"{_EXAMPLE_MAGMA_STEM}.genes.out")
        .read_text()
        .splitlines()
    )
    kept_out = [genes_out[0]] + [ln for ln in genes_out[1:] if ln.split()[1] == "22"]
    (out_dir / "chr22.genes.out").write_text("\n".join(kept_out) + "\n")
    genes_raw = (
        (example_magma_dir / f"{_EXAMPLE_MAGMA_STEM}.genes.raw")
        .read_text()
        .splitlines()
    )
    kept_raw = genes_raw[:2] + [ln for ln in genes_raw[2:] if ln.split()[1] == "22"]
    (out_dir / "chr22.genes.raw").write_text("\n".join(kept_raw) + "\n")
    return out_prefix


def test_pops_lowmem_matches_stock(tmp_path: Path):
    """The low-memory POPs reimplementation must produce PoPS scores and coefficients
    identical (to numerical precision) to stock pops.py. Runs both on the tarball's
    example munged features and a small chr22 MAGMA subset and compares outputs."""
    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run([POPS_SOURCE_EXTRACTED], incremental_save=True)
        source_asset = result[POPS_SOURCE_EXTRACTED.asset_id]
        assert isinstance(source_asset, DirectoryAsset)
        source_dir = source_asset.path

        magma_prefix = _write_chr22_magma_subset(
            source_dir / "example" / "data" / "magma_scores", tmp_path / "magma"
        )
        common_args = [
            "--gene_annot_path",
            str(source_dir / GENE_ANNOT_RELATIVE_PATH),
            "--feature_mat_prefix",
            str(source_dir / "example" / "data" / "features_munged" / "pops_features"),
            "--num_feature_chunks",
            "2",
            "--magma_prefix",
            str(magma_prefix),
            "--control_features_path",
            str(source_dir / CONTROL_FEATURES_RELATIVE_PATH),
        ]
        stock_prefix = tmp_path / "stock"
        lowmem_prefix = tmp_path / "lowmem"
        invoke_pops_script(
            source_dir,
            POPS_SCRIPT_NAME,
            common_args + ["--out_prefix", str(stock_prefix)],
        )
        invoke_pops_lowmem_script(
            source_dir, common_args + ["--out_prefix", str(lowmem_prefix)]
        )

        stock_preds = pd.read_csv(f"{stock_prefix}.preds", sep="\t").set_index("ENSGID")
        lowmem_preds = pd.read_csv(f"{lowmem_prefix}.preds", sep="\t").set_index(
            "ENSGID"
        )
        assert stock_preds.index.equals(lowmem_preds.index)
        max_pred_diff = np.abs(
            stock_preds["PoPS_Score"] - lowmem_preds["PoPS_Score"]
        ).max()
        assert max_pred_diff < 1e-9, (
            f"PoPS scores diverge (max abs diff {max_pred_diff})"
        )

        stock_coefs = pd.read_csv(f"{stock_prefix}.coefs", sep="\t").set_index(
            "parameter"
        )["beta"]
        lowmem_coefs = pd.read_csv(f"{lowmem_prefix}.coefs", sep="\t").set_index(
            "parameter"
        )["beta"]
        assert stock_coefs["SELECTED_CV_ALPHA"] == lowmem_coefs["SELECTED_CV_ALPHA"]
        meta = ["METHOD", "SELECTED_CV_ALPHA", "BEST_CV_SCORE"]
        stock_betas = stock_coefs[~stock_coefs.index.isin(meta)].astype(float)
        lowmem_betas = lowmem_coefs[~lowmem_coefs.index.isin(meta)].astype(float)
        assert set(stock_betas.index) == set(lowmem_betas.index)
        max_coef_diff = np.abs(
            stock_betas - lowmem_betas.reindex(stock_betas.index)
        ).max()
        assert max_coef_diff < 1e-9, (
            f"coefficients diverge (max abs diff {max_coef_diff})"
        )


def test_pops_source_utils_files_present(tmp_path: Path):
    """The gene-annotation and control-feature files POPs runs depend on must exist
    in the extracted source at the relative paths pops_utils references."""
    from mecfs_bio.build_system.task.pops.pops_utils import (
        CONTROL_FEATURES_RELATIVE_PATH,
        GENE_ANNOT_RELATIVE_PATH,
    )

    info_store = tmp_path / "info_store.yaml"
    asset_root = tmp_path / "asset_store"

    with log_on_error(info_store):
        asset_root.mkdir(parents=True, exist_ok=True)
        test_runner = SimpleRunner(
            tracer=ImoHasher.with_xxhash_128(),
            info_store=info_store,
            asset_root=asset_root,
        )
        result = test_runner.run([POPS_SOURCE_EXTRACTED], incremental_save=True)
        source_asset = result[POPS_SOURCE_EXTRACTED.asset_id]
        assert isinstance(source_asset, DirectoryAsset)
        assert (source_asset.path / "pops.py").is_file()
        assert (source_asset.path / "munge_feature_directory.py").is_file()
        assert (source_asset.path / GENE_ANNOT_RELATIVE_PATH).is_file()
        assert (source_asset.path / CONTROL_FEATURES_RELATIVE_PATH).is_file()
