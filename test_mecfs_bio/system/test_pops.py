"""
System test for the POPs tasks: munge features, then run POPs.

Uses the small example data bundled inside the POPs source tarball (the
PASS_Schizophrenia example from the POPs tutorial): 8 raw feature files over ~18k
genes and a matching MAGMA gene-analysis result. Exercises both
PopsMungeFeatureDirectoryTask and PopsRunTask end to end and verifies that a
non-empty .preds table is produced.
"""

from pathlib import Path

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
    POPS_OUTPUT_STEM_NAME,
    POPS_PREDS_SUFFIX,
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
