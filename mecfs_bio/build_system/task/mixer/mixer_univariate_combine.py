import json
import shutil
import tempfile
from pathlib import Path
from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.mixer_task import MIXER_FIT_JSON_PATTERN, MIXER_TEST_JSON_PATTERN
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer_figures
from mecfs_bio.build_system.wf.base_wf import WF


_CONTAINER_AGGREGATION_DIR = Path("/container_agg")

COMBINED_FIT_FILENAME_PREFIX = "trait1.fit"
COMBINED_TEST_FILENAME_PREFIX = "trait1.test"


@frozen
class MixerRunSource:
    task: Task
    rep:int

@frozen
class MixerUnivariateCombine(Task):
    mixer_source_runs: Sequence[MixerRunSource]
    _meta: Meta
    trait_name:str
    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [item.task for item in self.mixer_source_runs]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:

        with tempfile.TemporaryDirectory() as tmpdir_name:
            tmp_path = Path(tmpdir_name)
            agg_mounts = {tmp_path.resolve(): _CONTAINER_AGGREGATION_DIR}
            for source_run in self.mixer_source_runs:
                source_asset = fetch(source_run.task.asset_id)
                assert isinstance(source_asset, DirectoryAsset)
                shutil.copytree(source_asset.path, tmp_path, dirs_exist_ok=True)
                _edit_json_to_fix_trait_path(tmp_path/MIXER_FIT_JSON_PATTERN.replace("@", str(source_run.rep)), trait_name=self.trait_name )
                _edit_json_to_fix_trait_path(tmp_path/MIXER_TEST_JSON_PATTERN.replace("@", str(source_run.rep)), trait_name=self.trait_name )


            invoke_mixer_figures(
                args=[
                    "combine",
                    "--json",
                    str(_CONTAINER_AGGREGATION_DIR/MIXER_FIT_JSON_PATTERN),
                    "--out",
                    str(_CONTAINER_AGGREGATION_DIR / COMBINED_FIT_FILENAME_PREFIX),
                ],
                extra_mounts=agg_mounts,
            )

            invoke_mixer_figures(
                args=[
                    "combine",
                    "--json",
                    str(_CONTAINER_AGGREGATION_DIR/MIXER_TEST_JSON_PATTERN),
                    "--out",
                    str(_CONTAINER_AGGREGATION_DIR / COMBINED_TEST_FILENAME_PREFIX),
                ],
                extra_mounts=agg_mounts,
            )
            Path(tmp_path / (COMBINED_FIT_FILENAME_PREFIX + ".json")).rename(scratch_dir / (COMBINED_FIT_FILENAME_PREFIX + ".json"))
            Path(tmp_path / (COMBINED_TEST_FILENAME_PREFIX + ".json")).rename(scratch_dir / (COMBINED_TEST_FILENAME_PREFIX + ".json"))
            return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
            asset_id: str,
        mixer_source_runs: Sequence[MixerRunSource],
            trait_name:str
    ):
        assert len(mixer_source_runs)>= 1
        source_meta =mixer_source_runs[0].task.meta
        if isinstance(source_meta, SimpleDirectoryMeta):
            meta= SimpleDirectoryMeta(
                id=AssetId(asset_id),
            )
        elif isinstance(source_meta, ResultDirectoryMeta):
            meta =ResultDirectoryMeta(
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
            trait_name=trait_name,
        )



def _edit_json_to_fix_trait_path(json_path: Path, trait_name:str):
    with open(json_path, 'r') as f:
        data = json.load(f)

    # Overwrite the dynamic path with a single, static dummy path
    # this is needed, because MIXER's combine functionality checks that trait1_file is the same among inputs
    assert "trait1_file" in data["options"]
    data["options"]["trait1_file"] = f"standardized_path/{trait_name}"
    with open(json_path, 'w') as f:
        json.dump(data, f, indent=4)
