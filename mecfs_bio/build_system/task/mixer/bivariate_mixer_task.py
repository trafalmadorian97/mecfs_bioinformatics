import tempfile
from pathlib import Path
from typing import Callable, Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.mixer_task import MixerDataSource, PreformattedMixerDataSource, MixerTask, \
    CONTAINER_REF_DIR, prepare_mixer_trait_input_file, get_mixer_extract_args, MIXER_FIT_JSON_PATTERN
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer
from mecfs_bio.build_system.wf.base_wf import WF

_TRAIT_1_STATS = Path("/trait_1_stats")
_TRAIT_2_STATS = Path("/trait_2_stats")
_TRAIT_1_RESULT_DIR = Path("/trait_1_result")
_TRAIT_2_RESULT_DIR = Path("/trait_2_result")
_OUT_DIR = Path("/bivariate_mixer_output")




MIXER_BIVARIATE_FIT_PREFIX_PATTERN = "mixer_bivariate.fit.@"
MIXER_BIVARIATE_FIT_JSON_PATTERN = MIXER_BIVARIATE_FIT_PREFIX_PATTERN+".json"
MIXER_BIVARIATE_TEST_PREFIX_PATTERN = "mixer_bivariate.test.@"
@frozen
class BivariateMixerTask(Task):
    trait_1_source: MixerDataSource | PreformattedMixerDataSource
    trait_2_source: MixerDataSource | PreformattedMixerDataSource
    reference_data_directory_task: Task
    trait_1_univariate_task: MixerTask
    trait_2_univariate_task: MixerTask
    extract_file_pattern_gen: Callable[[int], str] | None
    extra_args: Sequence[str] = tuple()
    chr_to_use_arg: str | None = None
    ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
    bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
    threads: int = 4


    def __attrs_post_init__(self):
        assert len(self.trait_1_univariate_task.reps_to_perform)==1
        assert len(self.trait_2_univariate_task.reps_to_perform)==1
        assert self.trait_1_univariate_task.reps_to_perform == self.trait_2_univariate_task.reps_to_perform

    @property
    def rep(self)-> int:
        return self.trait_1_univariate_task.reps_to_perform[0]

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [
            self.trait_1_source.task,
            self.trait_2_source.task,
            self.trait_1_univariate_task,
            self.trait_2_univariate_task,
            self.reference_data_directory_task,
        ]
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        chr_args = (
            ["--chr2use", self.chr_to_use_arg]
            if self.chr_to_use_arg is not None
            else []
        )
        reference_dir_asset = fetch(self.reference_data_directory_task.asset_id)
        assert isinstance(reference_dir_asset, DirectoryAsset)
        with tempfile.TemporaryDirectory() as tempdir:
            temp_path =Path(tempdir)
            trait_1_stats_path=prepare_mixer_trait_input_file(
                source=self.trait_1_source,
                fetch=fetch,
                temp_dir=temp_path,
            )

            trait_2_stats_path=prepare_mixer_trait_input_file(
                source=self.trait_2_source,
                fetch=fetch,
                temp_dir=temp_path,
            )

            trait_1_asset = fetch(self.trait_1_univariate_task.asset_id)
            trait_2_asset = fetch(self.trait_2_univariate_task.asset_id)
            assert isinstance(trait_1_asset, DirectoryAsset)
            assert isinstance(trait_2_asset, DirectoryAsset)
            trait_1_result_dir = trait_1_asset.path
            trait_2_result_dir = trait_2_asset.path

            mounts = {reference_dir_asset.path.resolve(): CONTAINER_REF_DIR,
                      trait_1_stats_path.resolve():_TRAIT_1_STATS,
                      trait_2_stats_path.resolve():_TRAIT_2_STATS,
                      trait_1_result_dir.resolve():_TRAIT_1_RESULT_DIR,
                      trait_2_result_dir.resolve():_TRAIT_2_RESULT_DIR,
                      scratch_dir.resolve(): _OUT_DIR,
                      }

            common_args = [
                "--ld-file",
                str(CONTAINER_REF_DIR / self.ld_file_pattern),
                "--bim-file",
                str(CONTAINER_REF_DIR / self.bim_file_pattern),
                "--threads",
                str(self.threads),
            ]

            extract_args = get_mixer_extract_args(
                extract_file_pattern_gen=self.extract_file_pattern_gen,
                rep=self.rep,
                reference_dir_path=CONTAINER_REF_DIR,
            )
            invoke_mixer(
                ["fit2"]
                + common_args
                + extract_args
                +chr_args
                + list(self.extra_args)
                + [
                    "--trait1-file",
                    str(_TRAIT_1_STATS),
                    "--trait2-file",
                    str(_TRAIT_2_STATS),
                    "--trait1-params",
                    str(_TRAIT_1_RESULT_DIR/MIXER_FIT_JSON_PATTERN.replace("@",str(self.rep))),
                    "--trait2-params",
                    str(_TRAIT_2_RESULT_DIR/MIXER_FIT_JSON_PATTERN.replace("@",str(self.rep))),
                    # "--load-params",
                    # fit1_trait1_out_path_prefix + ".json",
                    "--out",
                    str(_OUT_DIR/MIXER_BIVARIATE_FIT_PREFIX_PATTERN.replace("@",str(self.rep))),
                    # test1_out_path_prefix,
                ],
                extra_mounts=mounts,
            )
            invoke_mixer(
                ["test2"]
                + common_args
                +chr_args
                +[
                    "--trait1-file",
                    str(_TRAIT_1_STATS),
                    "--trait2-file",
                    str(_TRAIT_2_STATS),
                    "--load-params",
                    str(_OUT_DIR/MIXER_BIVARIATE_FIT_JSON_PATTERN.replace("@",str(self.rep))),
                        "--out",
                    str(_OUT_DIR/MIXER_BIVARIATE_TEST_PREFIX_PATTERN.replace("@",str(self.rep))),
                ],
                extra_mounts=mounts,
            )
