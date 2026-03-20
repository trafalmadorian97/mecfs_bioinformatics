import os
import shutil
import tempfile
from pathlib import Path, PurePath
from typing import Callable, Sequence

from attrs import frozen
from structlog import get_logger

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.mixer.mixer_task import (
    CONTAINER_REF_DIR,
    MIXER_FIT_JSON_PATTERN,
    MixerDataSource,
    MixerTask,
    PreformattedMixerDataSource,
    get_mixer_extract_args,
    prepare_mixer_trait_input_file,
)
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer
from mecfs_bio.build_system.wf.base_wf import WF

logger = get_logger()

MIXER_BIVARIATE_FIT_PREFIX_PATTERN = "mixer_bivariate.fit.@"
MIXER_BIVARIATE_FIT_JSON_PATTERN = MIXER_BIVARIATE_FIT_PREFIX_PATTERN + ".json"
MIXER_BIVARIATE_TEST_PREFIX_PATTERN = "mixer_bivariate.test.@"
MIXER_BIVARIATE_TEST_JSON_PATTERN = MIXER_BIVARIATE_TEST_PREFIX_PATTERN + ".json"


@frozen
class BivariateMixerTask(Task):
    """
    Bivariate (cross-trait) MiXeR analysis.

    Runs fit2 and test2 steps, which require completed univariate fit1 results
    for both traits. Each BivariateMixerTask handles a single rep.

    See: O.Frei et al., Bivariate causal mixture model quantifies polygenic overlap
    between complex traits beyond genetic correlation, Nature Communications, 2019.

    Added with help from Claude
    """

    _meta: Meta
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
        assert len(self.trait_1_univariate_task.reps_to_perform) == 1
        assert len(self.trait_2_univariate_task.reps_to_perform) == 1
        assert (
            self.trait_1_univariate_task.reps_to_perform
            == self.trait_2_univariate_task.reps_to_perform
        )

    @property
    def rep(self) -> int:
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
        ref_mounts = {reference_dir_asset.path.resolve(): CONTAINER_REF_DIR}

        trait_1_asset = fetch(self.trait_1_univariate_task.asset_id)
        trait_2_asset = fetch(self.trait_2_univariate_task.asset_id)
        assert isinstance(trait_1_asset, DirectoryAsset)
        assert isinstance(trait_2_asset, DirectoryAsset)

        rep = self.rep
        fit_json = MIXER_FIT_JSON_PATTERN.replace("@", str(rep))
        fit_prefix = MIXER_BIVARIATE_FIT_PREFIX_PATTERN.replace("@", str(rep))
        test_prefix = MIXER_BIVARIATE_TEST_PREFIX_PATTERN.replace("@", str(rep))

        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tempdir:
            tmp_path = Path(tempdir).relative_to(os.getcwd())

            # Prepare trait sumstats files in temp dir
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

            # Copy univariate fit params into temp dir
            trait1_params_name = f"trait1_params.{rep}.json"
            trait2_params_name = f"trait2_params.{rep}.json"
            shutil.copy(
                str(trait_1_asset.path / fit_json),
                str(tmp_path / trait1_params_name),
            )
            shutil.copy(
                str(trait_2_asset.path / fit_json),
                str(tmp_path / trait2_params_name),
            )

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
                rep=rep,
                reference_dir_path=reference_dir_asset.path,
            )

            bivar_fit_out = str(tmp_path / fit_prefix)
            invoke_mixer(
                ["fit2"]
                + common_args
                + extract_args
                + chr_args
                + list(self.extra_args)
                + [
                    "--trait1-file",
                    str(trait_1_stats_path),
                    "--trait2-file",
                    str(trait_2_stats_path),
                    "--trait1-params",
                    str(tmp_path / trait1_params_name),
                    "--trait2-params",
                    str(tmp_path / trait2_params_name),
                    "--out",
                    bivar_fit_out,
                ],
                extra_mounts=ref_mounts,
            )

            bivar_test_out = str(tmp_path / test_prefix)
            invoke_mixer(
                ["test2"]
                + common_args
                + chr_args
                + [
                    "--trait1-file",
                    str(trait_1_stats_path),
                    "--trait2-file",
                    str(trait_2_stats_path),
                    "--load-params",
                    bivar_fit_out + ".json",
                    "--out",
                    bivar_test_out,
                ],
                extra_mounts=ref_mounts,
            )

            # Move outputs to scratch_dir
            for suffix in (".json", ".log"):
                fit_file = Path(bivar_fit_out + suffix)
                if fit_file.exists():
                    fit_file.rename(scratch_dir / fit_file.name)
                test_file = Path(bivar_test_out + suffix)
                if test_file.exists():
                    test_file.rename(scratch_dir / test_file.name)

            return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        trait_1_source: MixerDataSource | PreformattedMixerDataSource,
        trait_2_source: MixerDataSource | PreformattedMixerDataSource,
        ref_data_directory_task: Task,
        trait_1_univariate_task: MixerTask,
        trait_2_univariate_task: MixerTask,
        extract_file_pattern_gen: Callable[[int], str] | None,
        extra_args: Sequence[str] = tuple(),
        chr_args: str | None = None,
        ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld",
        bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim",
        threads: int = 4,
    ):
        source_meta = trait_1_source.task.meta
        meta: Meta
        if isinstance(source_meta, FilteredGWASDataMeta):
            meta = ResultDirectoryMeta(
                id=asset_id,
                trait="multi_trait",
                project="polygenic_overlap",
                sub_dir=PurePath("analysis") / "bivariate_mixer",
            )
        elif isinstance(source_meta, SimpleDirectoryMeta):
            meta = SimpleDirectoryMeta(
                id=AssetId(asset_id),
            )
        else:
            raise ValueError(f"Unknown meta {source_meta}")
        return cls(
            meta=meta,
            trait_1_source=trait_1_source,
            trait_2_source=trait_2_source,
            reference_data_directory_task=ref_data_directory_task,
            trait_1_univariate_task=trait_1_univariate_task,
            trait_2_univariate_task=trait_2_univariate_task,
            extract_file_pattern_gen=extract_file_pattern_gen,
            ld_file_pattern=ld_file_pattern,
            bim_file_pattern=bim_file_pattern,
            chr_to_use_arg=chr_args,
            extra_args=extra_args,
            threads=threads,
        )
