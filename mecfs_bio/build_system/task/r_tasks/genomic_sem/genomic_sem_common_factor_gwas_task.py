import gc
import tempfile
from pathlib import Path, PurePath
from typing import Sequence

from attrs import frozen
from rpy2 import robjects as ro
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
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_task import (
    GenomicSEMConfig,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_user_gwas_task import (
    COMMON_FACTOR_GWAS_FILENAME,
    GWAS_RESULTS_SUBDIR,
    GenomicSEMGWASRunConfig,
    GenomicSEMGWASSumstatsSource,
    GenomicSEMSumstatsConfig,
    _prepare_gwas_inputs,
    _r_to_pandas,
    _resolve_file_path,
    _resolve_ld_path,
    _validate_sources,
    logger,
)
from mecfs_bio.build_system.wf.base_wf import WF


@frozen
class GenomicSEMCommonFactorGWASTask(Task):
    """
    Run GenomicSEM common factor GWAS: munge → ldsc → sumstats →
    common factor GWAS. Output is a parquet of per-SNP common factor effects.
    """

    meta: Meta
    sources: Sequence[GenomicSEMGWASSumstatsSource]
    ld_ref_task: Task
    hapmap_snps_task: Task
    sumstats_ref_task: Task
    munge_config: GenomicSEMConfig = GenomicSEMConfig()
    sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig()
    run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig()

    def __attrs_post_init__(self):
        _validate_sources(self.sources)

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
            logger.debug("Running GenomicSEM::commonfactorGWAS")
            result = _run_common_factor_gwas(
                gsem=gsem, covstruc=covstruc, snps=snps, config=self.run_config
            )
            _save_common_factor_gwas_output(result, scratch_dir)

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
        munge_config: GenomicSEMConfig = GenomicSEMConfig(),
        sumstats_config: GenomicSEMSumstatsConfig = GenomicSEMSumstatsConfig(),
        run_config: GenomicSEMGWASRunConfig = GenomicSEMGWASRunConfig(),
    ) -> "GenomicSEMCommonFactorGWASTask":
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
            munge_config=munge_config,
            sumstats_config=sumstats_config,
            run_config=run_config,
        )


def _run_common_factor_gwas(*, gsem, covstruc, snps, config: GenomicSEMGWASRunConfig):
    """
    See https://github.com/GenomicSEM/GenomicSEM/wiki/4.-Common-Factor-GWAS
    """
    kwargs = dict(
        covstruc=covstruc,
        SNPs=snps,
        estimation=config.estimation,
        parallel=config.parallel,
        GC=config.gc_correction,
        toler=config.toler,
        SNPSE=config.snp_se,
        smooth_check=config.smooth_check,
    )
    if config.cores is not None:
        kwargs["cores"] = config.cores
    return gsem.commonfactorGWAS(**kwargs)


def _save_common_factor_gwas_output(result, scratch_dir: Path) -> Path:
    out_dir = scratch_dir / GWAS_RESULTS_SUBDIR
    out_dir.mkdir(parents=True, exist_ok=True)
    df = _r_to_pandas(result)
    out_path = out_dir / COMMON_FACTOR_GWAS_FILENAME
    df.to_parquet(out_path, index=False)
    logger.debug(f"Wrote common factor GWAS sumstats to {out_path} ({len(df)} rows)")
    return out_path
