"""
Task to combine bivariate MiXeR run outputs to produce a single result.
"""

import json
import shutil
import tempfile
from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.bivariate_mixer_task import BivariateMixerTask, MIXER_BIVARIATE_FIT_JSON_PATTERN, \
    MIXER_BIVARIATE_TEST_JSON_PATTERN
from mecfs_bio.build_system.task.mixer.mixer_task import (
    MIXER_FIT_JSON_PATTERN,
    MIXER_TEST_JSON_PATTERN,
)
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer_figures
from mecfs_bio.build_system.wf.base_wf import WF

_CONTAINER_AGGREGATION_DIR = Path("/container_agg")

BIVARIATE_COMBINED_FIT_FILENAME_PREFIX = "mixer_bivariate.fit"
BIVARIATE_COMBINED_TEST_FILENAME_PREFIX = "mixer_bivariate.test"


@frozen
class BivariateMixerRunSource:
    task: BivariateMixerTask
    rep: int


@frozen
class MixerBivariateCombine(Task):
    """
    Task to combine Bivariate MiXeR run outputs to produce a single result.

    The MiXeR authors have split up the genetic variants in their reference panel into 20 random subsets.
    The recommended MiXeR workflow is to run MiXeR on your GWAS data using each of these 20 random subsets, then combine the results.

    Use this task to combine the results of these runs
    """

    mixer_source_runs: Sequence[BivariateMixerRunSource]
    _meta: Meta


    @property
    def trait_1_name(self)->str:
        return self.mixer_source_runs[0].task.trait_1_source.alias

    @property
    def trait_2_name(self)->str:
        return self.mixer_source_runs[0].task.trait_2_source.alias

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.mixer_source_runs]

    def __attrs_post_init__(self):
        assert len(self.mixer_source_runs) >= 1
        for run in self.mixer_source_runs:
            assert run.task.trait_1_source.alias==self.mixer_source_runs[0].task.trait_1_source.alias
            assert run.task.trait_2_source.alias==self.mixer_source_runs[0].task.trait_2_source.alias


    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmp_path = Path(tmpdir_name)
            agg_mounts = {tmp_path.resolve(): _CONTAINER_AGGREGATION_DIR}
            for source_run in self.mixer_source_runs:
                source_asset = fetch(source_run.task.asset_id)
                assert isinstance(source_asset, DirectoryAsset)
                shutil.copytree(source_asset.path, tmp_path, dirs_exist_ok=True)
                # _edit_json_to_fix_trait_path(
                #     tmp_path / MIXER_FIT_JSON_PATTERN.replace("@", str(source_run.rep)),
                #     trait_1_name=self.trait_1_name,
                #     trait_2_name=self.trait_2_name,
                # )
                # _edit_json_to_fix_trait_path(
                #     tmp_path
                #     / MIXER_TEST_JSON_PATTERN.replace("@", str(source_run.rep)),
                #     trait_1_name=self.trait_1_name,
                #     trait_2_name=self.trait_2_name,
                # )

            invoke_mixer_figures(
                args=[
                    "combine",
                    "--json",
                    str(_CONTAINER_AGGREGATION_DIR /MIXER_BIVARIATE_FIT_JSON_PATTERN ),
                    "--out",
                    str(_CONTAINER_AGGREGATION_DIR / BIVARIATE_COMBINED_FIT_FILENAME_PREFIX),
                ],
                extra_mounts=agg_mounts,
            )

            invoke_mixer_figures(
                args=[
                    "combine",
                    "--json",
                    str(_CONTAINER_AGGREGATION_DIR / MIXER_BIVARIATE_TEST_JSON_PATTERN  ),
                    "--out",
                    str(_CONTAINER_AGGREGATION_DIR / BIVARIATE_COMBINED_TEST_FILENAME_PREFIX),
                ],
                extra_mounts=agg_mounts,
            )
            Path(tmp_path / (BIVARIATE_COMBINED_FIT_FILENAME_PREFIX + ".json")).rename(
                scratch_dir / (BIVARIATE_COMBINED_FIT_FILENAME_PREFIX + ".json")
            )
            Path(tmp_path / (BIVARIATE_COMBINED_TEST_FILENAME_PREFIX + ".json")).rename(
                scratch_dir / (BIVARIATE_COMBINED_TEST_FILENAME_PREFIX + ".json")
            )
            return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls, asset_id: str, mixer_source_runs: Sequence[BivariateMixerRunSource],
    ):
        assert len(mixer_source_runs) >= 1
        source_meta = mixer_source_runs[0].task.meta
        meta: Meta
        if isinstance(source_meta, SimpleDirectoryMeta):
            meta = SimpleDirectoryMeta(
                id=AssetId(asset_id),
            )
        elif isinstance(source_meta, ResultDirectoryMeta):
            meta = ResultDirectoryMeta(
                id=AssetId(asset_id),
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=source_meta.sub_dir,
            )
        else:
            raise NotImplementedError(f"Unknown source meta: {source_meta}")
        return cls(
            mixer_source_runs=mixer_source_runs,
            meta=meta,
        )


def _edit_json_to_fix_trait_path(json_path: Path, trait_1_name: str, trait_2_name: str):
    with open(json_path) as f:
        data = json.load(f)

    # Overwrite the dynamic path with a single, static dummy path
    # this is needed, because MIXER's combine functionality checks that trait1_file is the same among inputs
    assert "trait1_file" in data["options"]
    assert "trait2_file" in data["options"]
    data["options"]["trait1_file"] = f"standardized_path/{trait_1_name}"
    data["options"]["trait2_file"] = f"standardized_path/{trait_2_name}"
    with open(json_path, "w") as f:
        json.dump(data, f, indent=4)
