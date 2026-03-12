import os

from tqdm import tqdm
import tempfile
from typing import Sequence, Mapping, Callable

import narwhals
import polars as pl
from pathlib import Path, PurePath

import narwhals as nw
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import PhenotypeInfo, QuantPhenotype, \
    BinaryPhenotypeSampleInfo
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL, GWASLAB_CHROM_COL, GWASLAB_POS_COL, \
    GWASLAB_EFFECT_ALLELE_COL, GWASLAB_NON_EFFECT_ALLELE_COL, GWASLAB_BETA_COL, GWASLAB_SE_COL
from mecfs_bio.util.subproc.run_command import execute_command

from structlog import get_logger

MIXER_RSID_COL  ="RSID"
MIXER_CHROM_COL = "CHR"
MIXER_POS_COL = "POS"
MIXER_EFFECT_ALLELE_COL = "EffectAllele"
MIXER_NON_EFFECT_ALLELE_COL = "OtherAllele"
MIXER_EFFECTIVE_SAMPLE_SIZE = "N"
MIXER_Z_SCORE_COL = "Z"


MIXER_VERSION = "2.2.1"

logger = get_logger()


def _default_extract_file_pattern_gen(rep: int)->str:
    return "1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep{}.snps".format(rep)

@frozen
class MixerDataSource:
    """
    A source for data for use in Mixer
    """

    task: Task
    alias: str
    sample_info: PhenotypeInfo
    pipe: DataProcessingPipe = IdentityPipe()

    @property
    def asset_id(self) -> AssetId:
        return self.task.asset_id

@frozen
class BivariateMode:
    trait_2_source:MixerDataSource

@frozen
class UnivariateMode:
    pass

MixerMode = BivariateMode | UnivariateMode

@frozen
class MixerTask(Task):
    _meta: Meta
    trait_1_source: MixerDataSource
    # trait_2_source: MixerDataSource
    mixer_mode: MixerMode
    reference_data_directory_task: Task
    extract_file_pattern_gen:Callable[[int],str] | None
    ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
    bim_file_pattern:str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
    chr_to_use_arg:str|None=None
    threads: int=4
    num_reps: int=20


    def __attrs_post_init__(self):
        _invoke_mixer("--version",{})

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        result= [
            self.trait_1_source.task,self.reference_data_directory_task
        ]
        if isinstance(self.mixer_mode, BivariateMode):
            result.append(self.mixer_mode.trait_2_source.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        reference_dir_asset = fetch(self.reference_data_directory_task.asset_id)
        assert isinstance(reference_dir_asset, DirectoryAsset)
        with tempfile.TemporaryDirectory(dir=os.getcwd()) as tmpdir:
            tmp_path = Path(tmpdir).relative_to(os.getcwd())
            trait1_path = _prep_summary_statistics_for_mixer(
                sumstats_dataframe_task=self.trait_1_source.task,
                fetch=fetch,
                pipe=self.trait_1_source.pipe,
                phenotype=self.trait_1_source.sample_info,
                name=self.trait_1_source.alias,
                temp_dir=tmp_path,
            )
            assert trait1_path.is_file()
            if isinstance(self.mixer_mode, BivariateMode):
                trait2_path = _prep_summary_statistics_for_mixer(
                    sumstats_dataframe_task=self.mixer_mode.trait_2_source.task,
                    fetch=fetch,
                    pipe=self.mixer_mode.trait_2_source.pipe,
                    phenotype=self.mixer_mode.trait_2_source.sample_info,
                    name=self.mixer_mode.trait_2_source.alias,
                    temp_dir=tmp_path,
                )
                assert trait2_path.is_file()
            else:
                trait2_path = None
            common_args =  ["--ld-file", str(reference_dir_asset.path/self.ld_file_pattern), "--bim-file",
                            str(reference_dir_asset.path/self.bim_file_pattern),
                            "--threads", str(self.threads)]

            for rep in tqdm(range(1, self.num_reps + 1)):

                if self.extract_file_pattern_gen is not None:
                    extract_file = reference_dir_asset.path / self.extract_file_pattern_gen(rep)
                    assert extract_file.is_file()
                    extract_args = ["--extract", str(extract_file)]
                else:
                    extract_args=[]
                if self.chr_to_use_arg is not None:
                    chr_args=["--chr2use", self.chr_to_use_arg]
                else:
                    chr_args=[]
                fit1_trait1_out_path_prefix = str(tmp_path/f"trait1.fit.{rep}")
                _invoke_mixer(
                  ["fit1"]+  common_args + extract_args+ chr_args+ ["--trait1-file", str(trait1_path),
                                                 "--out",
                                                 str(fit1_trait1_out_path_prefix)
                                                 ],
                    extra_mounts={
                    }
                )
                fit1_trait1_out_json = Path(fit1_trait1_out_path_prefix + ".json")
                fit1_trait1_out_log = Path(fit1_trait1_out_path_prefix + ".log")
                fit1_trait1_out_json.rename(scratch_dir/f"trait1.fit.{rep}.json")
                fit1_trait1_out_log.rename(scratch_dir/f"trait1.fit.{rep}.log")

                # to add later: handling of bivariate case



            return DirectoryAsset(scratch_dir)
    @classmethod
    def create(cls,
               asset_id:str,
               trait_1_source: MixerDataSource,
               trait_2_source: MixerDataSource,

               ce_data_directory_task: Task,
        ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld",

    bim_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim",
    extract_file_pattern_gen: Callable[[int],str]|None = _default_extract_file_pattern_gen ,
    threads: int = 4,
    num_reps: int = 20,
               chr_args: str|None=None
    ):
        meta =  ResultDirectoryMeta(
            id=asset_id,
            trait = "multi_trait",
            project="polygenic_overlap",
            sub_dir=PurePath("mixer")
        )
        return cls(
            meta=meta,
            trait_1_source=trait_1_source,
            trait_2_source=trait_2_source,
            reference_data_directory_task=ce_data_directory_task,
            ld_file_pattern=ld_file_pattern,
            bim_file_pattern=bim_file_pattern,
            extract_file_pattern_gen=extract_file_pattern_gen,
            threads=threads,
            num_reps=num_reps,
            chr_to_use_arg=chr_args,
        )


def _prep_summary_statistics_for_mixer(
    sumstats_dataframe_task:    Task,
    fetch: Fetch,
    temp_dir: Path,
    pipe: DataProcessingPipe,
    phenotype: PhenotypeInfo,
    name:str) -> Path:
    asset = fetch(sumstats_dataframe_task.asset_id)
    frame= pipe.process(scan_dataframe_asset(asset, meta=sumstats_dataframe_task.meta))
    frame = frame.with_columns(
        nw.col(GWASLAB_RSID_COL).alias(MIXER_RSID_COL),
        nw.col(GWASLAB_CHROM_COL).alias(MIXER_CHROM_COL),
        nw.col(GWASLAB_POS_COL).alias(MIXER_POS_COL),
        nw.col( GWASLAB_EFFECT_ALLELE_COL).alias(MIXER_EFFECT_ALLELE_COL),
        nw.col(GWASLAB_NON_EFFECT_ALLELE_COL).alias(MIXER_NON_EFFECT_ALLELE_COL),
        (nw.col(GWASLAB_BETA_COL)/   nw.col(GWASLAB_SE_COL)  ).alias(MIXER_Z_SCORE_COL)
    ).select(
        MIXER_RSID_COL,
        MIXER_CHROM_COL,
        MIXER_POS_COL,
        MIXER_EFFECT_ALLELE_COL,
        MIXER_NON_EFFECT_ALLELE_COL,
        MIXER_Z_SCORE_COL
    )
    if isinstance(phenotype, QuantPhenotype):
        assert phenotype.total_sample_size is not None
        frame = frame.with_columns(
            narwhals.lit(phenotype.total_sample_size).alias(MIXER_EFFECTIVE_SAMPLE_SIZE),
        )
    elif isinstance(phenotype, BinaryPhenotypeSampleInfo):
        frame = frame.with_columns(
            narwhals.lit(
                phenotype.effective_sample_size
            ).alias(MIXER_EFFECTIVE_SAMPLE_SIZE)
        )
    out_path = temp_dir/name
    frame.collect().to_pandas().to_csv(out_path, index=False, sep="\t")
    return out_path


SETUP_MIXER_DOCKER = [
    'export',f'MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:{MIXER_VERSION} python /tools/mixer/precimed/mixer.py";'
]

def _get_docker_command(extra_mounts:Mapping[Path, Path])->list[str]:
    inner_docker_command = "docker run --shm-size=2g -v $PWD:/home"
    for key, value in extra_mounts.items():
        inner_docker_command+= f" -v {str(key)}:{str(value)}"
    inner_docker_command+=" -w /home"
    return[
        f'export',f'DOCKER_RUN="{inner_docker_command}";'
    ]


def _invoke_mixer(
        args: Sequence[str]|str,
        extra_mounts:Mapping[Path, Path],
):
    if isinstance(args, str):
        args = [args]
    execute_command(
       _get_docker_command(extra_mounts=extra_mounts)+SETUP_MIXER_DOCKER+["${MIXER_PY}"] +list(args)
    )

