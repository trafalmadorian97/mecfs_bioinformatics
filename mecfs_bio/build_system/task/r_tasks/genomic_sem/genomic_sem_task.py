"""
Run a basic GenomicSEM workflow: munge, multivariable LDSC, and a user-specified
factor model expressed in lavaan syntax.

See: https://github.com/GenomicSEM/GenomicSEM/wiki/3.-Genome%E2%80%90wide-Models

GenomicSEM is an R library, so it is accessed through Python via rpy2.

Shared configuration/dataclasses live in `_genomic_sem_config`, rpy2-free input
helpers in `_genomic_sem_inputs`, and the R-calling helpers in
`_genomic_sem_r_bridge`. This module holds only the `GenomicSEMTask` class.
"""

import gc
import tempfile
from pathlib import Path, PurePath
from typing import Sequence

import rpy2.robjects as ro
import structlog
from attrs import frozen
from rpy2.robjects.packages import importr

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.result_directory_meta import ResultDirectoryMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    MULTI_TRAIT,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    LAVAAN_MODEL_FILENAME,
    LDSC_LOG_PREFIX,
    MUNGED_SUBDIR,
    GenomicSEMConfig,
    GenomicSEMSumstatsSource,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_inputs import (
    _get_prevs,
    _get_sample_size,
    _write_munge_input,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_r_bridge import (
    _run_ldsc,
    _run_munge,
    _save_ldsc_outputs,
    _save_model_outputs,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class GenomicSEMTask(Task):
    """
    Run the basic GenomicSEM workflow:

    1. Munge the supplied summary statistics into the format expected by LDSC.
    2. Run multivariable LDSC on the munged sumstats to obtain a genetic
       covariance matrix.
    3. Fit a factor model expressed in the lavaan language to the LDSC output.

    Inputs:
      - ld_ref_task: produces a directory containing reference LD scores in standard form.
      - hapmap_snps_task: produces a file listing the HapMap3 SNPs to keep during munging.
      - sources: list of tasks producing tabular GWAS summary statistics in
        gwaslab format.
      - factor_model: a factor model description in lavaan syntax.

    Output: a directory asset containing munged sumstats, LDSC outputs, and
    factor model results plus logs from intermediate steps.
    """

    meta: Meta
    sources: Sequence[GenomicSEMSumstatsSource]
    ld_ref_task: Task
    hapmap_snps_task: Task
    factor_model: str
    config: GenomicSEMConfig = GenomicSEMConfig()

    def __attrs_post_init__(self):
        assert len(self.sources) >= 2, "GenomicSEM requires at least two traits"
        aliases = [s.alias for s in self.sources]
        assert len(set(aliases)) == len(aliases), (
            f"Trait aliases must be unique. Got: {aliases}"
        )

    @property
    def deps(self) -> list["Task"]:
        result: list[Task] = [self.ld_ref_task, self.hapmap_snps_task]
        for source in self.sources:
            result.append(source.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gsem = importr("GenomicSEM")

        ld_asset = fetch(self.ld_ref_task.asset_id)
        assert isinstance(ld_asset, DirectoryAsset)
        # Absolute path so values remain correct after R chdirs during munge.
        ld_path = str(ld_asset.path.resolve())

        hapmap_asset = fetch(self.hapmap_snps_task.asset_id)
        assert isinstance(hapmap_asset, FileAsset)
        hapmap_path = str(hapmap_asset.path.resolve())

        munged_dir = scratch_dir / MUNGED_SUBDIR
        munged_dir.mkdir(parents=True, exist_ok=True)

        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            input_files = []
            trait_names = []
            sample_sizes: list[float] = []
            sample_prevs: list[float] = []
            population_prevs: list[float] = []
            for source in self.sources:
                input_path = _write_munge_input(
                    source=source, fetch=fetch, tmp_dir=tmp_dir
                )
                input_files.append(str(input_path))
                trait_names.append(source.alias)
                sample_sizes.append(_get_sample_size(source))
                samp_prev, pop_prev = _get_prevs(source.sample_info)
                sample_prevs.append(samp_prev)
                population_prevs.append(pop_prev)

            logger.debug(f"Running GenomicSEM::munge on {len(input_files)} files")
            _run_munge(
                gsem=gsem,
                input_files=input_files,
                hapmap_path=hapmap_path,
                trait_names=trait_names,
                sample_sizes=sample_sizes,
                output_dir=munged_dir,
                info_filter=self.config.info_filter,
                maf_filter=self.config.maf_filter,
            )

            munged_paths = [
                str(munged_dir / f"{name}.sumstats.gz") for name in trait_names
            ]
            for path in munged_paths:
                assert Path(path).is_file(), f"Munged sumstats not found at {path}"

            logger.debug("Running GenomicSEM::ldsc")
            covstruc = _run_ldsc(
                gsem=gsem,
                munged_paths=munged_paths,
                trait_names=trait_names,
                sample_prevs=sample_prevs,
                population_prevs=population_prevs,
                ld_path=ld_path,
                ld_file_basename_prefix=self.config.ld_file_filename_pattern,
                ldsc_log_prefix=str(scratch_dir / LDSC_LOG_PREFIX),
            )

            _save_ldsc_outputs(covstruc=covstruc, scratch_dir=scratch_dir)

            (scratch_dir / LAVAAN_MODEL_FILENAME).write_text(self.factor_model)

            logger.debug("Running GenomicSEM::usermodel")
            model_result = gsem.usermodel(
                covstruc=covstruc,
                estimation=self.config.estimation,
                model=self.factor_model,
                CFIcalc=self.config.cfi_calc,
                std_lv=self.config.std_lv,
                fix_resid=self.config.fix_resid,
            )
            _save_model_outputs(model_result=model_result, scratch_dir=scratch_dir)

        gc.collect()
        ro.r("gc()")
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sources: Sequence[GenomicSEMSumstatsSource],
        ld_ref_task: Task,
        hapmap_snps_task: Task,
        factor_model: str,
        config: GenomicSEMConfig = GenomicSEMConfig(),
    ) -> "GenomicSEMTask":
        meta = ResultDirectoryMeta(
            id=AssetId(asset_id),
            trait=MULTI_TRAIT,
            project="genomic_sem",
            sub_dir=PurePath("analysis"),
        )
        return cls(
            meta=meta,
            sources=sources,
            ld_ref_task=ld_ref_task,
            hapmap_snps_task=hapmap_snps_task,
            factor_model=factor_model,
            config=config,
        )
