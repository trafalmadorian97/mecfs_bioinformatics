"""
Unit tests for the .create() meta derivation of the SBayesRC and polypwas tasks.

These exercise the trait/project derivation and output meta wiring without running
Docker.
"""

from pathlib import PurePath

import pytest

from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.meta.simple_file_meta import SimpleFileMeta
from mecfs_bio.build_system.task.fake_task import FakeTask
from mecfs_bio.build_system.task.polypwas.polypwas_assoc_task import PolypwasAssocTask
from mecfs_bio.build_system.task.polypwas.polypwas_train_task import (
    POLYPWAS_WEIGHTS_EXTENSION,
    PolypwasTrainTask,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_data_source import (
    PreformattedSBayesRCDataSource,
)
from mecfs_bio.build_system.task.sbayesrc.sbayesrc_task import SBayesRCTask


def _gwas_summary_source(trait: str) -> PreformattedSBayesRCDataSource:
    return PreformattedSBayesRCDataSource(
        task=FakeTask(
            meta=GWASSummaryDataFileMeta(
                id=f"{trait}_ma",
                trait=trait,
                project="proj",
                sub_dir="processed",
                project_path=None,
                extension=".ma",
            )
        ),
        filename=None,
        alias=trait,
    )


def _ld_task() -> FakeTask:
    return FakeTask(meta=SimpleDirectoryMeta("ld_ref"))


def test_sbayesrc_create_derives_result_directory_meta():
    task = SBayesRCTask.create(
        asset_id="sbrc_run",
        gwas_source=_gwas_summary_source("ldl"),
        ld_reference_directory_task=_ld_task(),
    )
    assert isinstance(task.meta, ResultDirectoryMeta)
    assert task.meta.trait == "ldl"
    assert task.meta.project == "proj"
    assert task.meta.sub_dir == PurePath("analysis") / "sbayesrc"


def test_polypwas_train_create_derives_weights_meta():
    task = PolypwasTrainTask.create(
        asset_id="train_run",
        pqtl_source=_gwas_summary_source("angptl3"),
        ld_reference_directory_task=_ld_task(),
    )
    assert isinstance(task.meta, FilteredGWASDataMeta)
    assert task.meta.trait == "angptl3"
    assert task.meta.extension == POLYPWAS_WEIGHTS_EXTENSION


def test_polypwas_assoc_create_derives_result_table_meta():
    task = PolypwasAssocTask.create(
        asset_id="assoc_run",
        weights_task=FakeTask(meta=SimpleFileMeta("weights")),
        gwas_source=_gwas_summary_source("ldl"),
        ld_reference_directory_task=_ld_task(),
        gene_info_task=FakeTask(meta=SimpleFileMeta("gene_info")),
    )
    assert isinstance(task.meta, ResultTableMeta)
    assert task.meta.trait == "ldl"
    assert task.meta.extension == ".tsv"
    assert task.meta.read_spec is not None


def _assoc_task(weights_task, weights_filename):
    return PolypwasAssocTask.create(
        asset_id="assoc_run",
        weights_task=weights_task,
        gwas_source=_gwas_summary_source("ldl"),
        ld_reference_directory_task=_ld_task(),
        gene_info_task=FakeTask(meta=SimpleFileMeta("gene_info")),
        weights_filename=weights_filename,
    )


def test_assoc_weights_filename_invariant():
    file_weights = FakeTask(meta=SimpleFileMeta("weights"))
    dir_weights = FakeTask(meta=SimpleDirectoryMeta("weights_dir"))

    # Valid combinations construct fine.
    _assoc_task(file_weights, None)
    _assoc_task(dir_weights, "weights.wgts.gz")

    with pytest.raises(AssertionError):
        _assoc_task(file_weights, "weights.wgts.gz")
    with pytest.raises(AssertionError):
        _assoc_task(dir_weights, None)
