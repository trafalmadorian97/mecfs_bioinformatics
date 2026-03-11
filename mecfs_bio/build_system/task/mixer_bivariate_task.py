import tempfile
from typing import Sequence

import narwhals
import polars as pl
from pathlib import Path

import narwhals as nw
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
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

MIXER_RSID_COL  ="RSID"
MIXER_CHROM_COL = "CHR"
MIXER_POS_COL = "POS"
MIXER_EFFECT_ALLELE_COL = "EffectAllele"
MIXER_NON_EFFECT_ALLELE_COL = "OtherAllele"
MIXER_EFFECTIVE_SAMPLE_SIZE = "N"
MIXER_Z_SCORE_COL = "Z"


MIXER_VERSION = "2.2.1"

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
class MixerBivariateTask(Task):
    _meta: Meta
    trait_1_source: MixerDataSource
    trait_2_source: MixerDataSource
    reference_data_directory_task: Task
    ld_file_pattern: str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.run4.ld"
    bim_file_pattern:str = "1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim"
    rep_file_pattern:str = r"1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.{}.snps"
    threads: int=4
    num_reps: int=20


    def __attrs_post_init__(self):
        _invoke_mixer("--version")

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [
            self.trait_1_source.task,self.trait_2_source.task,self.reference_data_directory_task
        ]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            trait1_path = _prep_summary_statistics_for_mixer(
                sumstats_dataframe_task=self.trait_1_source.task,
                fetch=fetch,
                pipe=self.trait_1_source.pipe,
                phenotype=self.trait_1_source.sample_info,
                name=self.trait_1_source.alias,
                temp_dir=tmp_path,
            )
            trait2_path = _prep_summary_statistics_for_mixer(
                sumstats_dataframe_task=self.trait_2_source.task,
                fetch=fetch,
                pipe=self.trait_2_source.pipe,
                phenotype=self.trait_2_source.sample_info,
                name=self.trait_2_source.alias,
                temp_dir=tmp_path,
            )
            common_args =  ["--ld-file", self.ld_file_pattern, "--bim-file", self.bim_file_pattern, "--threads", str(self.threads)]

            for rep in range(1, self.num_reps + 1):
                extract_args = ["--extract",self.rep_file_pattern.format(rep)]
                _invoke_mixer(
                    common_args + extract_args+ ["--trait1-file", trait1_path,
                                                 "--out",
                                                 f"trait1.fit.{rep}"
                                                 ],
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
    frame.collect().to_pandas().to_csv(out_path, index=False)
    return out_path


SETUP_DOCKER_COMMAND = [
    'export','DOCKER_RUN="docker run -v $PWD:/home -w /home";'
]
SETUP_MIXER_DOCKER = [
    'export',f'MIXER_PY="$DOCKER_RUN ghcr.io/precimed/gsa-mixer:{MIXER_VERSION} python /tools/mixer/precimed/mixer.py";'
]



def _invoke_mixer(
        args: Sequence[str]|str,
):
    if isinstance(args, str):
        args = [args]
    execute_command(
       SETUP_DOCKER_COMMAND+SETUP_MIXER_DOCKER+["${MIXER_PY}"] +list(args)
    )

