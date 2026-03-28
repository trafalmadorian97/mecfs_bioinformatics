from pathlib import Path, PurePath

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.mixer_bivariate_combine import (
    BIVARIATE_COMBINED_FIT_FILENAME_PREFIX,
    BIVARIATE_COMBINED_TEST_FILENAME_PREFIX,
    MixerBivariateCombine,
)
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer_figures
from mecfs_bio.build_system.wf.base_wf import WF

_CONTAINER_PLOT_DIR = Path("/container_plot")
_CONTAINER_COMBINED_DIR = Path("/container_combine")

BIVARIATE_OUTPUT_PREFIX = "bivariate_results"


@frozen
class MixerBivariateSummarizeResultsTask(Task):
    """
    Task to take the results of a Bivarirate MiXeR run and produce tables and plots.
    """

    _meta: Meta
    combine_task: MixerBivariateCombine
    override_trait_1_name: str | None = None
    override_trait_2_name: str | None = None

    @property
    def trait_1_name(self) -> str:
        return self.combine_task.trait_1_name

    @property
    def trait_2_name(self) -> str:
        return self.combine_task.trait_2_name

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.combine_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        source_asset = fetch(self.combine_task.asset_id)
        assert isinstance(source_asset, DirectoryAsset)

        plots_mounts = {
            scratch_dir.resolve(): _CONTAINER_PLOT_DIR,
            source_asset.path.resolve(): _CONTAINER_COMBINED_DIR,
        }
        trait_1_name = self.override_trait_1_name or self.trait_1_name
        trait_2_name = self.override_trait_2_name or self.trait_2_name
        invoke_mixer_figures(
            args=[
                "two",
                "--json-test",
                str(
                    _CONTAINER_COMBINED_DIR
                    / (BIVARIATE_COMBINED_TEST_FILENAME_PREFIX + ".json")
                ),
                "--json-fit",
                str(
                    _CONTAINER_COMBINED_DIR
                    / (BIVARIATE_COMBINED_FIT_FILENAME_PREFIX + ".json")
                ),
                "--out",
                str(_CONTAINER_PLOT_DIR / (BIVARIATE_OUTPUT_PREFIX)),
                "--statistic",
                "mean std",
                "--ext",
                "png",
                "--trait1",
                _normalize_name(trait_1_name),
                "--trait2",
                _normalize_name(trait_2_name),
            ],
            extra_mounts=plots_mounts,
        )
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        combine_task: MixerBivariateCombine,
        override_trait_1_name: str | None = None,
        override_trait_2_name: str | None = None,
    ):
        source_meta = combine_task.meta
        meta: Meta
        if isinstance(source_meta, ResultDirectoryMeta):
            meta = ResultDirectoryMeta(
                trait=source_meta.trait,
                project=source_meta.project,
                id=AssetId(asset_id),
                sub_dir=PurePath("analysis") / "mixer_results",
            )
        elif isinstance(source_meta, SimpleDirectoryMeta):
            meta = SimpleDirectoryMeta(
                id=AssetId(asset_id),
            )
        else:
            raise ValueError(f"Unknown meta {source_meta}")
        return cls(
            meta=meta,
            combine_task=combine_task,
            override_trait_1_name=override_trait_1_name,
            override_trait_2_name=override_trait_2_name,
        )


def _normalize_name(name: str) -> str:
    return name.replace("_", "-")
