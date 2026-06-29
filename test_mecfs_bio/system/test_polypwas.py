"""
System tests for the polypwas tasks (train + assoc).

These are slow-tier tests: they run polypwas download-example, which fetches the
HapMap3 LD reference (~3 GB) plus the bundled ANGPTL3 + LDL example, and the train
test additionally drives SBayesRC through the Docker image.  They are intended for
scheduled / on-demand runs, not the default fast suite.

Reference (polypwas README, Option A/B): the bundled ANGPTL3 + LDL example yields
CIS_Z ~ 17.21 and TRANS_Z ~ -45.68 when using the published pre-trained weights.
"""

import subprocess
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.rebuilder.verifying_trace_rebuilder.tracer.imohash import (
    ImoHasher,
)
from mecfs_bio.build_system.runner.simple_runner import SimpleRunner
from mecfs_bio.build_system.task.polypwas.polypwas_assoc_task import PolypwasAssocTask
from mecfs_bio.build_system.task.polypwas.polypwas_train_task import PolypwasTrainTask
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    PreformattedSBayesRCDataSource,
)
from test_mecfs_bio.system.util import (
    StampedExternalDirectoryTask,
    StampedExternalFileTask,
    log_on_error,
)

_PROJECT = "polypwas_example"


def _download_example(tmp_path: Path, include_weights: bool) -> Path:
    """Run polypwas download-example in tmp_path and return its data directory."""
    args = ["polypwas", "download-example"]
    if include_weights:
        args.append("--include-weights")
    subprocess.run(args, cwd=tmp_path, check=True)
    return tmp_path / "data"


def _ld_task(data_dir: Path) -> StampedExternalDirectoryTask:
    return StampedExternalDirectoryTask(
        meta=SimpleDirectoryMeta(id="polypwas_example_ld"),
        external_path=data_dir / "ldm" / "ukbEUR_HM3",
    )


def _ma_source(
    data_dir: Path, filename: str, trait: str
) -> PreformattedSBayesRCDataSource:
    task = StampedExternalFileTask(
        meta=GWASSummaryDataFileMeta(
            id=f"polypwas_example_{trait}_ma",
            trait=trait,
            project=_PROJECT,
            sub_dir="processed",
            project_path=None,
            extension=".ma.gz",
        ),
        external_path=data_dir / "examples" / filename,
    )
    return PreformattedSBayesRCDataSource(task=task, filename=filename, alias=trait)


def _make_runner(tmp_path: Path) -> SimpleRunner:
    asset_root = tmp_path / "asset_store"
    asset_root.mkdir(parents=True, exist_ok=True)
    return SimpleRunner(
        tracer=ImoHasher.with_xxhash_128(),
        info_store=tmp_path / "info_store.yaml",
        asset_root=asset_root,
    )


def test_polypwas_assoc_pretrained(tmp_path: Path):
    """
    Deterministic assoc check using the published pre-trained ANGPTL3 weights.

    polypwas assoc needs no R, so this exercises PolypwasAssocTask end-to-end
    against the readme's reference Z-scores.
    """
    data_dir = _download_example(tmp_path, include_weights=True)

    weights_task = StampedExternalFileTask(
        meta=SimpleFileMeta(id="polypwas_example_angptl3_weights"),
        external_path=data_dir / "examples" / "angptl3.wgts.gz",
    )
    gene_info_task = StampedExternalFileTask(
        meta=SimpleFileMeta(id="polypwas_example_angptl3_gene_info"),
        external_path=data_dir / "examples" / "angptl3.gene.tsv",
    )
    assoc_task = PolypwasAssocTask.create(
        asset_id="polypwas_assoc_angptl3_ldl_pretrained",
        weights_task=weights_task,
        gwas_source=_ma_source(data_dir, "ldl.ma.gz", "ldl"),
        ld_reference_directory_task=_ld_task(data_dir),
        gene_info_task=gene_info_task,
    )

    info_store = tmp_path / "info_store.yaml"
    with log_on_error(info_store):
        result = _make_runner(tmp_path).run([assoc_task], incremental_save=True)
        assert result is not None
        out_path = result[assoc_task.asset_id].path
        table = pl.read_csv(out_path, separator="\t")
        assert {"ID", "CIS_Z", "TRANS_Z"}.issubset(set(table.columns))
        row = table.filter(pl.col("ID") == "ANGPTL3").row(0, named=True)
        assert abs(row["CIS_Z"] - 17.21) < 1.0
        assert abs(row["TRANS_Z"] - (-45.68)) < 2.0


def test_polypwas_train_then_assoc(tmp_path: Path):
    """
    Train ANGPTL3 weights via the SBayesRC Docker path, then run assoc.

    SBayesRC's random seed makes the exact Z-scores vary, so we assert the cis
    signal is strongly positive and the trans signal strongly negative rather than
    pinning exact values.
    """
    data_dir = _download_example(tmp_path, include_weights=False)

    train_task = PolypwasTrainTask.create(
        asset_id="polypwas_train_angptl3",
        pqtl_source=_ma_source(data_dir, "angptl3.ma.gz", "angptl3"),
        ld_reference_directory_task=_ld_task(data_dir),
        threads=4,
    )
    gene_info_task = StampedExternalFileTask(
        meta=SimpleFileMeta(id="polypwas_example_angptl3_gene_info"),
        external_path=data_dir / "examples" / "angptl3.gene.tsv",
    )
    assoc_task = PolypwasAssocTask.create(
        asset_id="polypwas_assoc_angptl3_ldl_trained",
        weights_task=train_task,
        gwas_source=_ma_source(data_dir, "ldl.ma.gz", "ldl"),
        ld_reference_directory_task=_ld_task(data_dir),
        gene_info_task=gene_info_task,
    )

    info_store = tmp_path / "info_store.yaml"
    with log_on_error(info_store):
        result = _make_runner(tmp_path).run(
            [train_task, assoc_task], incremental_save=True
        )
        assert result is not None
        out_path = result[assoc_task.asset_id].path
        table = pl.read_csv(out_path, separator="\t")
        row = table.filter(pl.col("ID") == "ANGPTL3").row(0, named=True)
        assert row["CIS_Z"] > 5.0
        assert row["TRANS_Z"] < -10.0
