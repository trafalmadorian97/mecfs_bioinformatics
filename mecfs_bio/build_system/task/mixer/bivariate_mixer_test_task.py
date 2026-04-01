import os
import shutil
import tempfile
from pathlib import Path, PurePath
from typing import Callable, Sequence

from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.bivariate_mixer_task import (
    MIXER_BIVARIATE_FIT_PREFIX_PATTERN,
    MIXER_BIVARIATE_TEST_PREFIX_PATTERN,
    BivariateMixerTask,
)
from mecfs_bio.build_system.task.mixer.mixer_task import (
    CONTAINER_REF_DIR,
    MIXER_FIT_JSON_PATTERN,
    MixerDataSource,
    PreformattedMixerDataSource,
    get_mixer_extract_args,
    prepare_mixer_trait_input_file,
)
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class BivariateMixerTestTask(Task):
    """
    Task to run bivariate MiXeR test evaluation on a model produced by bivariate MiXeR fit
    """

    meta: Meta
    trait_1_source: MixerDataSource | PreformattedMixerDataSource
    trait_2_source: MixerDataSource | PreformattedMixerDataSource
    reference_data_directory_task: Task
    fit_task: BivariateMixerTask
    extra_args: Sequence[str] = tuple()
    chr_to_use_arg: str | None = None
    ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
    bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
    threads: int = 4
    extract_file_pattern_gen: Callable[[int], str] | None = None

    @property
    def rep(self) -> int:
        return self.fit_task.rep

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [
            self.trait_1_source.task,
            self.trait_2_source.task,
            self.reference_data_directory_task,
            self.fit_task,
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
        ref_mounts = {reference_dir_asset.path.resolve(): CONTAINER_REF_DIR}

        rep = self.rep
        fit_json = MIXER_FIT_JSON_PATTERN.replace("@", str(rep))
        bivariate_fit_prefix = MIXER_BIVARIATE_FIT_PREFIX_PATTERN.replace("@", str(rep))
        bivariate_fit_json = bivariate_fit_prefix + ".json"
        bivariate_test_prefix = MIXER_BIVARIATE_TEST_PREFIX_PATTERN.replace(
            "@", str(rep)
        )

        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            tmp_path = Path(tempdir).relative_to(os.getcwd())

            fit_asset = fetch(self.fit_task.asset_id)
            assert isinstance(fit_asset, DirectoryAsset)
            bivariate_fit_json_path = tmp_path / bivariate_fit_json
            shutil.copy(fit_asset.path / (bivariate_fit_json), bivariate_fit_json_path)

            trait_1_stats_path = prepare_mixer_trait_input_file(
                source=self.trait_1_source,
                fetch=fetch,
                temp_dir=tmp_path,
            )
            trait_2_stats_path = prepare_mixer_trait_input_file(
                source=self.trait_2_source,
                fetch=fetch,
                temp_dir=tmp_path,
            )
            common_args = [
                "--ld-file",
                str(CONTAINER_REF_DIR / self.ld_file_pattern),
                "--bim-file",
                str(CONTAINER_REF_DIR / self.bim_file_pattern),
                "--threads",
                str(self.threads),
            ]

            bivar_test_out_str = str(tmp_path / bivariate_test_prefix)

            extra_test_args = []
            if self.extract_file_pattern_gen is not None:
                extract_args = get_mixer_extract_args(
                    extract_file_pattern_gen=self.extract_file_pattern_gen,
                    rep=rep,
                    reference_dir_path=reference_dir_asset.path,
                )
                extra_test_args.extend(extract_args)
            invoke_mixer(
                ["test2"]
                + common_args
                + chr_args
                + extra_test_args
                + [
                    "--trait1-file",
                    str(trait_1_stats_path),
                    "--trait2-file",
                    str(trait_2_stats_path),
                    "--load-params",
                    str(bivariate_fit_json_path),
                    "--out",
                    bivar_test_out_str,
                ],
                extra_mounts=ref_mounts,
            )
            test_out_json_path = Path(bivar_test_out_str + ".json")
            test_out_log_path = Path(bivar_test_out_str + ".log")
            assert test_out_json_path.exists()
            assert test_out_log_path.exists()
            test_out_json_path.rename(tmp_path / test_out_json_path.name)
            test_out_log_path.rename(tmp_path / test_out_log_path.name)
            return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        trait_1_source: MixerDataSource | PreformattedMixerDataSource,
        trait_2_source: MixerDataSource | PreformattedMixerDataSource,
        ref_data_directory_task: Task,
        fit_task: BivariateMixerTask,
        extra_args: Sequence[str] = tuple(),
        chr_to_use_arg: str | None = None,
        ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld",
        bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim",
        threads: int = 4,
        extract_file_pattern_gen: Callable[[int], str] | None = None,
    ):
        meta = ResultDirectoryMeta(
            id=asset_id,
            trait="multi_trait",
            project="polygenic_overlap",
            sub_dir=PurePath("analysis") / "bivariate_mixer",
        )
        return cls(
            meta=meta,
            fit_task=fit_task,
            trait_1_source=trait_1_source,
            trait_2_source=trait_2_source,
            reference_data_directory_task=ref_data_directory_task,
            extract_file_pattern_gen=extract_file_pattern_gen,
            extra_args=extra_args,
            chr_to_use_arg=chr_to_use_arg,
            ld_file_pattern=ld_file_pattern,
            bim_file_pattern=bim_file_pattern,
            threads=threads,
        )
