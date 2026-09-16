"""
Task to read in a dataframe from disk and process it using the GWASLAB pipeline.
"""

import attrs
import structlog
from loguru import logger

from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.reference_meta.reference_file_meta import (
    ReferenceFileMeta,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.constants.genomic_coordinate_constants import GenomeBuild

logger = structlog.get_logger()
from pathlib import Path
from typing import Literal, Sequence

import gwaslab as gl
import narwhals
import pandas as pd
from attrs import frozen
from gwaslab.util.util_in_filter_value import _exclude_sexchr

from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.filtered_gwas_data_meta import FilteredGWASDataMeta
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.read_spec.read_dataframe import scan_dataframe_asset
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_STATUS_COL,
    GwaslabKnownFormat,
)

_VARIANT_KEY_COLUMNS = [
    GWASLAB_CHROM_COL,
    GWASLAB_POS_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
]

GenomeBuildMode = Literal["infer", "19", "38"]


@frozen
class GWASLabColumnSpecifiers:
    rsid: str | None = None
    snpid: str | None = None
    chrom: str | None = None
    pos: str | None = None
    ea: str | None = None
    nea: str | None = None
    OR: str | None = None
    se: str | None = None
    p: str | None = None
    info: str | None = None
    eaf: str | None = None
    neaf: str | None = None
    beta: str | None = None
    ncase: str | None = None
    ncontrol: str | None = None
    neff: str | None = None
    n: str | None = None
    or_95l: str | None = None
    or_95u: str | None = None
    chi_sq: str | None = None
    mlog10p: str | None = None
    maf: str | None = None

    def get_selection_pipe(self) -> SelectColPipe:
        fields = attrs.asdict(self)
        cols = []
        for v in fields.values():
            if v is not None:
                cols.append(v)
        return SelectColPipe(cols)


ValidGwaslabFormat = GwaslabKnownFormat | GWASLabColumnSpecifiers


def _validate_eaf_in_range(sumstats: gl.Sumstats) -> None:
    if GWASLAB_EFFECT_ALLELE_FREQ_COL not in sumstats.data.columns:
        return
    eaf = sumstats.data[GWASLAB_EFFECT_ALLELE_FREQ_COL]
    eaf_min = float(eaf.min(skipna=True))
    eaf_max = float(eaf.max(skipna=True))
    if pd.isna(eaf_min) or pd.isna(eaf_max):
        logger.debug("EAF column present, but entirely NAN")
        return
    assert 0 <= eaf_min <= 1 and 0 <= eaf_max <= 1, (
        f"Effect allele frequency column {GWASLAB_EFFECT_ALLELE_FREQ_COL!r} has values "
        f"outside the [0, 1] fraction range (observed min={eaf_min}, max={eaf_max}). "
    )


def _get_sumstats(
    x: narwhals.LazyFrame,
    fmt: ValidGwaslabFormat,
    drop_cols: Sequence[str],
) -> gl.Sumstats:
    if isinstance(fmt, GWASLabColumnSpecifiers):
        x = fmt.get_selection_pipe().process(x)
    x = x.drop(drop_cols)
    logger.debug("Collecting Narwhals Lazyframe and converting to pandas")
    collected_df = x.collect().to_pandas()
    if isinstance(fmt, GWASLabColumnSpecifiers):
        return gl.Sumstats(
            collected_df,
            rsid=fmt.rsid,
            snpid=fmt.snpid,
            chrom=fmt.chrom,
            pos=fmt.pos,
            ea=fmt.ea,
            nea=fmt.nea,
            OR=fmt.OR,
            se=fmt.se,
            p=fmt.p,
            info=fmt.info,
            eaf=fmt.eaf,
            neaf=fmt.neaf,
            beta=fmt.beta,
            ncase=fmt.ncase,
            ncontrol=fmt.ncontrol,
            neff=fmt.neff,
            n=fmt.n,
            OR_95L=fmt.or_95l,
            OR_95U=fmt.or_95u,
            chisq=fmt.chi_sq,
            mlog10p=fmt.mlog10p,
            maf=fmt.maf,
        )

    return gl.Sumstats(
        collected_df,
        fmt=fmt,
    )


@frozen
class GWASLabCreateSumstatsTask(Task):
    """
    Task that processes a DataFrame of GWAS summary statistics using the GWASLab pipeline.
    see: https://cloufield.github.io/gwaslab/SumstatsObject/

    """

    df_source_task: Task
    target_asset_id: AssetId
    basic_check: bool
    genome_build: GenomeBuildMode
    filter_hapmap3: bool = False
    filter_indels: bool = False
    filter_palindromic: bool = False
    exclude_hla: bool = False
    exclude_sexchr: bool = False
    liftover_to: GenomeBuild | None = None
    fmt: GwaslabKnownFormat | GWASLabColumnSpecifiers = "regenie"
    drop_col_list: Sequence[str] = tuple()
    pre_pipe: DataProcessingPipe = IdentityPipe()

    def __attrs_post_init__(self):
        assert self._source_meta is not None

    @property
    def _source_id(self) -> AssetId:
        return self.df_source_task.meta.asset_id

    @property
    def _source_meta(
        self,
    ) -> FilteredGWASDataMeta | GWASSummaryDataFileMeta | ReferenceFileMeta:
        meta = self.df_source_task.meta
        assert isinstance(
            meta, (FilteredGWASDataMeta, GWASSummaryDataFileMeta, ReferenceFileMeta)
        )
        return meta

    @property
    def meta(self) -> Meta:
        if isinstance(self._source_meta, ReferenceFileMeta):
            return GWASLabSumStatsMeta(
                id=self.target_asset_id,
                trait="reference_data_gwas",
                project=self._source_meta.group,
            )
        return GWASLabSumStatsMeta(
            id=self.target_asset_id,
            trait=self._source_meta.trait,
            project=self._source_meta.project,
            sub_dir="gwaslab_sumstats",
        )

    @property
    def deps(self) -> list["Task"]:
        return [self.df_source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        df = scan_dataframe_asset(asset=fetch(self._source_id), meta=self._source_meta)
        logger.debug("Applying pre-pipe")
        df = self.pre_pipe.process(df)
        logger.debug("Fetching source dataframe asset...")
        sumstats = _get_sumstats(df, self.fmt, drop_cols=self.drop_col_list)
        _validate_eaf_in_range(sumstats)
        transform_spec = GwasLabTransformSpec(
            basic_check=self.basic_check,
            genome_build=self.genome_build,
            filter_hapmap3=self.filter_hapmap3,
            filter_indels=self.filter_indels,
            filter_palindromic=self.filter_palindromic,
            exclude_hla=self.exclude_hla,
            exclude_sexchr=self.exclude_sexchr,
            liftover_to=self.liftover_to,
        )
        sumstats = transform_gwaslab_sumstats(sumstats, spec=transform_spec)
        out_path = scratch_dir / "pickled_sumstats.pickle"
        gl.dump_pickle(sumstats, path=str(out_path))
        return FileAsset(out_path)


def _sumstats_raise_on_error(sumstats: gl.Sumstats):
    error_status = sumstats.data[GWASLAB_STATUS_COL].astype(str) == "9999999"
    if error_status.any():
        raise ValueError("GWASLAB Error")
    if len(sumstats.data) == 0:
        raise ValueError("No rows survive GWASLAB quality control!")


@frozen
class GwasLabTransformSpec:
    basic_check: bool = True
    genome_build: GenomeBuildMode = "infer"
    filter_hapmap3: bool = False
    filter_indels: bool = False
    filter_palindromic: bool = False
    exclude_hla: bool = False
    exclude_sexchr: bool = False
    liftover_to: GenomeBuild | None = None


def _drop_duplicate_variant_keys(sumstats: gl.Sumstats) -> None:
    """
    Drop every row whose (CHR, POS, EA, NEA) key is shared with another row, in place.

    Liftover can map two distinct source variants to the same target coordinate with the
    same alleles but different statistics, leaving duplicate variant keys that are ambiguous
    for any downstream allele-keyed join. Because we cannot tell which source variant is
    correct at the collided coordinate, all colliding rows are dropped, not deduplicated to
    one. Genome-reference harmonization asserts key uniqueness, so these must be removed here.
    """
    key = [col for col in _VARIANT_KEY_COLUMNS if col in sumstats.data.columns]
    if len(key) < len(_VARIANT_KEY_COLUMNS):
        return
    duplicated = sumstats.data.duplicated(subset=key, keep=False)
    count = int(duplicated.sum())
    if count:
        logger.warning(
            "dropping rows with duplicate (CHR, POS, EA, NEA) keys", count=count
        )
        sumstats.data = sumstats.data[~duplicated].reset_index(drop=True)


def transform_gwaslab_sumstats(
    sumstats: gl.Sumstats,
    spec: GwasLabTransformSpec,
) -> gl.Sumstats:
    logger.debug("Running gwas summary statistics through gwaslab pipelines...")
    logger.debug(f"Initial sumstats has shape {sumstats.data.shape}")
    if spec.basic_check:
        sumstats.basic_check()
    if spec.genome_build == "infer":
        sumstats.fix_chr()
        sumstats.infer_build()
        build = sumstats.meta["gwaslab"]["genome_build"]
        forced_build = None
        print(f"Build is {build}")
    else:
        build = spec.genome_build
        forced_build = build

    if spec.liftover_to is not None and (build != spec.liftover_to):
        sumstats.liftover(to_build=spec.liftover_to, from_build=forced_build)
        sumstats.infer_build()
        assert sumstats.build == spec.liftover_to
    if spec.filter_hapmap3:
        sumstats.filter_hapmap3(inplace=True, build=forced_build)
    if spec.filter_indels:
        sumstats.filter_indel(inplace=True, mode="out")
    if spec.filter_palindromic:
        sumstats.filter_palindromic(inplace=True, mode="out")
    if spec.exclude_hla:
        sumstats.exclude_hla(inplace=True)
    if spec.exclude_sexchr:
        sumstats = _exclude_sexchr(sumstats)
    _drop_duplicate_variant_keys(sumstats)
    _sumstats_raise_on_error(sumstats)
    logger.debug(f"Finished gwaslab pipe.  Data has shape {sumstats.data.shape}")
    return sumstats
