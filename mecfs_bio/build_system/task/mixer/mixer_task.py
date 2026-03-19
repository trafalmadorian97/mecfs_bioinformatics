"""
Core task for fitting the MiXeR Gaussian mixture model to GWAS data

"""

import os
import shutil
import tempfile
from pathlib import Path, PurePath
from typing import Callable, Sequence

import narwhals as nw
from attrs import frozen
from structlog import get_logger
from tqdm import tqdm

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.meta.simple_directory_meta import SimpleDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    BinaryPhenotypeSampleInfo,
    PhenotypeInfo,
    QuantPhenotype,
)
from mecfs_bio.build_system.task.mixer.mixer_utils import invoke_mixer
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
    GWASLAB_SE_COL,
)

MIXER_RSID_COL = "RSID"
MIXER_CHROM_COL = "CHR"
MIXER_POS_COL = "POS"
MIXER_EFFECT_ALLELE_COL = "EffectAllele"
MIXER_NON_EFFECT_ALLELE_COL = "OtherAllele"
MIXER_EFFECTIVE_SAMPLE_SIZE = "N"
MIXER_Z_SCORE_COL = "Z"

logger = get_logger()


def _default_extract_file_pattern_gen(rep: int) -> str:
    return (
        f"1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep{rep}.snps"
    )


def _get_fit_filename_prefix(rep: int) -> str:
    return f"trait1.fit.{rep}"


MIXER_FIT_JSON_PATTERN = "trait1.fit.@.json"
MIXER_TEST_JSON_PATTERN = "trait1.test.@.json"


@frozen
class MixerDataSource:
    """
    A source for data for use in Mixer.
    The task should provide a dataframe in gwaslab format, which will be
    converted to MiXeR format (column renaming + Z = BETA/SE computation).
    """

    task: Task
    alias: str
    sample_info: PhenotypeInfo
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class PreformattedMixerDataSource:
    """
    A source for data that is already in MiXeR sumstats format
    (RSID, CHR, POS, EffectAllele, OtherAllele, Z, N).
    No gwaslab-to-mixer column conversion is performed.
    The task should provide a DirectoryAsset or FileAsset containing the sumstats file.
    """

    task: Task
    filename: str
    alias: str

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id


@frozen
class BivariateMode:
    trait_2_source: MixerDataSource


@frozen
class UnivariateMode:
    pass


MixerMode = BivariateMode | UnivariateMode

CONTAINER_REF_DIR = Path("/ref_data")


@frozen
class MixerTask(Task):
    """
    Core task to fit the MiXeR Gaussian mixture model to GWAS data

    See:
    Holland, Dominic, et al. "Beyond SNP heritability: Polygenicity and discoverability of phenotypes
    estimated with a univariate Gaussian mixture model." PLoS Genetics 16.5 (2020): e1008612.

    The MiXeR software is distributed via Docker image.  Before running MixerTask, verify that you have installed Docker
    and added yourself to the Docker user group.

    The MiXeR authors have split up the genetic variants in their reference panel into 20 random subsets.
    The recommended MiXeR workflow is to run MiXeR on your GWAS data using each of these 20 random subsets, then combine the results.
    Specify which of these random subsets to run using the reps_to_perform attribute.

    """

    _meta: Meta
    trait_1_source: MixerDataSource | PreformattedMixerDataSource
    reference_data_directory_task: Task
    extract_file_pattern_gen: Callable[[int], str] | None
    extra_args: Sequence[str] = tuple()
    ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
    bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
    chr_to_use_arg: str | None = None
    threads: int = 4
    reps_to_perform: Sequence[int] = tuple(range(1, 21))

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [
            self.trait_1_source.task,
            self.reference_data_directory_task,
        ]
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        reference_dir_asset = fetch(self.reference_data_directory_task.asset_id)
        assert isinstance(reference_dir_asset, DirectoryAsset)
        ref_mounts = {reference_dir_asset.path.resolve(): CONTAINER_REF_DIR}
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            tmp_path = Path(tmpdir).relative_to(os.getcwd())
            trait1_path = prepare_mixer_trait_input_file(
                source=self.trait_1_source,
                fetch=fetch,
                temp_dir=tmp_path,
            )
            assert trait1_path.is_file()

            common_args = [
                "--ld-file",
                str(CONTAINER_REF_DIR / self.ld_file_pattern),
                "--bim-file",
                str(CONTAINER_REF_DIR / self.bim_file_pattern),
                "--threads",
                str(self.threads),
            ]

            for rep in tqdm(self.reps_to_perform):
                extract_args = get_mixer_extract_args(
                    extract_file_pattern_gen=self.extract_file_pattern_gen,
                    rep=rep,
                    reference_dir_path=reference_dir_asset.path,
                )
                chr_args = (
                    ["--chr2use", self.chr_to_use_arg]
                    if self.chr_to_use_arg is not None
                    else []
                )
                fit1_trait1_out_path_prefix = str(tmp_path / f"trait1.fit.{rep}")
                invoke_mixer(
                    ["fit1"]
                    + common_args
                    + extract_args
                    + chr_args
                    + list(self.extra_args)
                    + [
                        "--trait1-file",
                        str(trait1_path),
                        "--out",
                        str(fit1_trait1_out_path_prefix),
                    ],
                    extra_mounts=ref_mounts,
                )

                test1_out_path_prefix = str(tmp_path / f"trait1.test.{rep}")
                invoke_mixer(
                    ["test1"]
                    + common_args
                    + extract_args
                    + chr_args
                    + [
                        "--trait1-file",
                        str(trait1_path),
                        "--load-params",
                        fit1_trait1_out_path_prefix + ".json",
                        "--out",
                        test1_out_path_prefix,
                    ],
                    extra_mounts=ref_mounts,
                )
                Path(test1_out_path_prefix + ".json").rename(
                    scratch_dir / f"trait1.test.{rep}.json"
                )
                Path(test1_out_path_prefix + ".log").rename(
                    scratch_dir / f"trait1.test.{rep}.log"
                )

                Path(fit1_trait1_out_path_prefix + ".json").rename(
                    scratch_dir / f"{_get_fit_filename_prefix(rep)}.json"
                )
                Path(fit1_trait1_out_path_prefix + ".log").rename(
                    scratch_dir / f"{_get_fit_filename_prefix(rep)}.log"
                )

            return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        trait_1_source: MixerDataSource | PreformattedMixerDataSource,
        ref_data_directory_task: Task,
        extra_args: Sequence[str] = tuple(),
        ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld",
        bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim",
        extract_file_pattern_gen: Callable[[int], str]
        | None = _default_extract_file_pattern_gen,
        threads: int = 4,
        reps_to_perform: Sequence[int] = tuple(range(1, 21)),
        chr_args: str | None = None,
    ):
        source_meta = trait_1_source.task.meta
        meta: Meta
        if isinstance(source_meta, FilteredGWASDataMeta):
            meta = ResultDirectoryMeta(
                id=asset_id,
                trait=source_meta.trait,
                project=source_meta.project,
                sub_dir=PurePath("analysis") / "mixer",
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
            reference_data_directory_task=ref_data_directory_task,
            ld_file_pattern=ld_file_pattern,
            bim_file_pattern=bim_file_pattern,
            extract_file_pattern_gen=extract_file_pattern_gen,
            threads=threads,
            reps_to_perform=reps_to_perform,
            chr_to_use_arg=chr_args,
            extra_args=extra_args,
        )


def prepare_mixer_trait_input_file(
    source: MixerDataSource | PreformattedMixerDataSource,
    fetch: Fetch,
    temp_dir: Path,
) -> Path:
    """Prepare a trait sumstats file in the temp dir, ready for MiXeR."""
    if isinstance(source, PreformattedMixerDataSource):
        source_asset = fetch(source.task.asset_id)
        if isinstance(source_asset, DirectoryAsset):
            source_file = source_asset.path / source.filename
        elif isinstance(source_asset, FileAsset):
            source_file = source_asset.path
        else:
            raise ValueError(f"Unexpected asset type: {type(source_asset)}")
        assert source_file.is_file(), f"Source file not found: {source_file}"
        dest = temp_dir / source.filename
        shutil.copy(str(source_file), str(dest))
        return dest
    elif isinstance(source, MixerDataSource):
        return _prep_summary_statistics_for_mixer(
            sumstats_dataframe_task=source.task,
            fetch=fetch,
            pipe=source.pipe,
            phenotype=source.sample_info,
            name=source.alias,
            temp_dir=temp_dir,
        )
    else:
        raise ValueError(f"Unexpected source type: {type(source)}")


def _prep_summary_statistics_for_mixer(
    sumstats_dataframe_task: Task,
    fetch: Fetch,
    temp_dir: Path,
    pipe: DataProcessingPipe,
    phenotype: PhenotypeInfo,
    name: str,
) -> Path:
    asset = fetch(sumstats_dataframe_task.asset_id)
    frame = pipe.process(scan_dataframe_asset(asset, meta=sumstats_dataframe_task.meta))
    frame = frame.with_columns(
        nw.col(GWASLAB_RSID_COL).alias(MIXER_RSID_COL),
        nw.col(GWASLAB_CHROM_COL).alias(MIXER_CHROM_COL),
        nw.col(GWASLAB_POS_COL).alias(MIXER_POS_COL),
        nw.col(GWASLAB_EFFECT_ALLELE_COL).alias(MIXER_EFFECT_ALLELE_COL),
        nw.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(MIXER_NON_EFFECT_ALLELE_COL),
        (nw.col(GWASLAB_BETA_COL) / nw.col(GWASLAB_SE_COL)).alias(MIXER_Z_SCORE_COL),
    ).select(
        MIXER_RSID_COL,
        MIXER_CHROM_COL,
        MIXER_POS_COL,
        MIXER_EFFECT_ALLELE_COL,
        MIXER_NON_EFFECT_ALLELE_COL,
        MIXER_Z_SCORE_COL,
    )
    if isinstance(phenotype, QuantPhenotype):
        assert phenotype.total_sample_size is not None
        frame = frame.with_columns(
            nw.lit(phenotype.total_sample_size).alias(MIXER_EFFECTIVE_SAMPLE_SIZE),
        )
    elif isinstance(phenotype, BinaryPhenotypeSampleInfo):
        frame = frame.with_columns(
            nw.lit(phenotype.effective_sample_size).alias(MIXER_EFFECTIVE_SAMPLE_SIZE)
        )
    out_path = temp_dir / name
    frame.collect().to_pandas().to_csv(out_path, index=False, sep="\t")
    return out_path


@frozen
class MixerLDGenerationTask(Task):
    """
    Implemented by Claude to facilitate testing.
    Generates .ld files from PLINK .bed/.bim/.fam files using mixer ld command.
    Copies all source files plus generated .ld files to the output directory.
    """

    _meta: Meta
    plink_data_task: Task
    chromosomes: tuple[int, ...]
    bfile_prefix_pattern: str = "g1000_eur_hm3_chr{chr}"
    r2min: str = "0.05"
    ldscore_r2min: str = "0.01"
    ld_window_kb: str = "10000"

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [self.plink_data_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        src_asset = fetch(self.plink_data_task.asset_id)
        assert isinstance(src_asset, DirectoryAsset)

        # Copy all source files to scratch_dir
        for f in src_asset.path.iterdir():
            if f.is_file():
                shutil.copy2(str(f), str(scratch_dir / f.name))

        # Generate .ld files using mixer ld, mounting source dir in Docker
        src_mounts = {src_asset.path.resolve(): CONTAINER_REF_DIR}
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            tmp_path = Path(tmpdir).relative_to(os.getcwd())
            for chri in self.chromosomes:
                bfile_prefix = self.bfile_prefix_pattern.format(chr=chri)
                ld_out = str(tmp_path / f"{bfile_prefix}.ld")
                invoke_mixer(
                    [
                        "ld",
                        "--bfile",
                        str(CONTAINER_REF_DIR / bfile_prefix),
                        "--r2min",
                        self.r2min,
                        "--ldscore-r2min",
                        self.ldscore_r2min,
                        "--out",
                        ld_out,
                        "--ld-window-kb",
                        self.ld_window_kb,
                    ],
                    extra_mounts=src_mounts,
                )
                ld_file = Path(ld_out)
                assert ld_file.is_file(), f"Expected LD file not generated: {ld_file}"
                shutil.copy2(str(ld_file), str(scratch_dir / ld_file.name))

        return DirectoryAsset(scratch_dir)


def get_mixer_extract_args(
    extract_file_pattern_gen: Callable[[int], str] | None,
    rep: int,
    reference_dir_path: Path,
) -> list[str]:
    if extract_file_pattern_gen is not None:
        extract_file = reference_dir_path / extract_file_pattern_gen(rep)
        assert extract_file.is_file()
        extract_args = [
            "--extract",
            str(CONTAINER_REF_DIR / extract_file_pattern_gen(rep)),
        ]
        return extract_args
    return []
