"""
GenomicSEM multivariate (user) GWAS: adds SNP effects to a user-supplied
lavaan model and produces per-SNP summary statistics for the latent factors.

Workflow:
  1. munge() — produces .sumstats.gz for LDSC.
  2. ldsc() — produces covstruc.
  3. sumstats() — aligns per-trait stats to a reference panel.
  4. userGWAS(covstruc, SNPs, model).

Reference: https://github.com/GenomicSEM/GenomicSEM/wiki/5.-Multivariate-GWAS

Shared configuration/dataclasses live in `_genomic_sem_config`, rpy2-free input
helpers in `_genomic_sem_inputs`, and the R-calling helpers in
`_genomic_sem_r_bridge`. This module holds only the `GenomicSEMUserGWASTask`
class.
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
    GenomicSEMConfig,
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsConfig,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_inputs import (
    _resolve_file_path,
    _resolve_ld_path,
    _sanitize_component_name,
    _validate_sources,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_r_bridge import (
    _prepare_gwas_inputs,
    _run_user_gwas,
    _save_user_gwas_outputs,
)
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class GenomicSEMUserGWASTask(Task):
    """
    Run GenomicSEM multivariate (user) GWAS: munge → ldsc → sumstats →
    userGWAS with a user-supplied lavaan model containing SNP effects.

    sub_components selects which model components are extracted per SNP and
    written to disk, e.g. ["F1~SNP", "F2~SNP"]. One parquet per component is
    written to gwas_results/.
    """

    meta: Meta
    sources: Sequence[GenomicSEMGWASSumstatsSource]
    ld_ref_task: Task
    hapmap_snps_task: Task
    sumstats_ref_task: Task
    factor_model: str
    sub_components: Sequence[str]
    munge_config: GenomicSEMConfig = GenomicSEMConfig()
    sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig()
    run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig()
    fix_measurement: bool = True
    q_snp: bool = False
    std_lv: bool = False

    def __attrs_post_init__(self):
        _validate_sources(self.sources)
        assert len(self.sub_components) >= 1, (
            "userGWAS requires at least one sub_components entry"
        )
        sanitized = [_sanitize_component_name(c) for c in self.sub_components]
        assert len(set(sanitized)) == len(sanitized), (
            f"sub_components map to colliding filenames: {sanitized}"
        )

    @property
    def deps(self) -> list[Task]:
        result: list[Task] = [
            self.ld_ref_task,
            self.hapmap_snps_task,
            self.sumstats_ref_task,
        ]
        for source in self.sources:
            result.append(source.task)
        return result

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        gsem = importr("GenomicSEM")
        ld_path = _resolve_ld_path(self.ld_ref_task, fetch, self.munge_config)
        hapmap_path = _resolve_file_path(self.hapmap_snps_task, fetch)
        sumstats_ref_path = _resolve_file_path(self.sumstats_ref_task, fetch)

        with tempfile.TemporaryDirectory() as tmp_dir_str:
            tmp_dir = Path(tmp_dir_str)
            covstruc, snps = _prepare_gwas_inputs(
                gsem=gsem,
                sources=self.sources,
                ld_path=ld_path,
                hapmap_path=hapmap_path,
                sumstats_ref_path=sumstats_ref_path,
                munge_config=self.munge_config,
                sumstats_config=self.sumstats_config,
                fetch=fetch,
                scratch_dir=scratch_dir,
                tmp_dir=tmp_dir,
            )
            (scratch_dir / LAVAAN_MODEL_FILENAME).write_text(self.factor_model)
            logger.debug("Running GenomicSEM::userGWAS")
            result = _run_user_gwas(
                gsem=gsem,
                covstruc=covstruc,
                snps=snps,
                model=self.factor_model,
                sub_components=self.sub_components,
                fix_measurement=self.fix_measurement,
                q_snp=self.q_snp,
                std_lv=self.std_lv,
                config=self.run_config,
            )
            _save_user_gwas_outputs(
                result=result,
                sub_components=self.sub_components,
                scratch_dir=scratch_dir,
            )

        gc.collect()
        ro.r("gc()")
        return DirectoryAsset(scratch_dir)

    @classmethod
    def create(
        cls,
        asset_id: str,
        sources: Sequence[GenomicSEMGWASSumstatsSource],
        ld_ref_task: Task,
        hapmap_snps_task: Task,
        sumstats_ref_task: Task,
        factor_model: str,
        sub_components: Sequence[str],
        munge_config: GenomicSEMConfig = GenomicSEMConfig(),
        sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig(),
        run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig(),
        fix_measurement: bool = True,
        q_snp: bool = False,
        std_lv: bool = False,
    ) -> "GenomicSEMUserGWASTask":
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
            sumstats_ref_task=sumstats_ref_task,
            factor_model=factor_model,
            sub_components=sub_components,
            munge_config=munge_config,
            sumstats_config=sumstats_config,
            run_config=run_config,
            fix_measurement=fix_measurement,
            q_snp=q_snp,
            std_lv=std_lv,
        )
