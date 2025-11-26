import structlog
from loguru import logger

from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe

logger = structlog.get_logger()
from pathlib import Path
from typing import Literal, Sequence

import gwaslab
import gwaslab as gl
import narwhals
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
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import (
    GWASLAB_STATUS_COL,
    GwaslabKnownFormat,
)
from mecfs_bio.build_system.wf.base_wf import WF

GenomeBuildMode = Literal["infer", "19", "38"]
GenomeBuild = Literal["19", "38"]


@frozen
class GWASLabVCFRef:
    name: str
    ref_alt_freq: str


@frozen
class HarmonizationOptions:
    """
    Options for the call to GWASLab's harmonize function
    """

    ref_infer: GWASLabVCFRef
    ref_seq: str
    cores: int
    check_ref_files: bool
    drop_missing_from_ref: bool


def _do_harmonization(
    sumstats: gl.Sumstats, basic_check: bool, options: HarmonizationOptions
):
    if options.check_ref_files:
        gwaslab.download_ref(name=options.ref_infer.name, overwrite=False)
        gwaslab.download_ref(name=options.ref_seq, overwrite=False)
    sumstats.harmonize(
        basic_check=basic_check,
        n_cores=options.cores,
        ref_seq=gl.get_path(options.ref_seq),
        ref_infer=gl.get_path(options.ref_infer.name),
        ref_alt_freq=options.ref_infer.ref_alt_freq,
    )
    if options.drop_missing_from_ref:
        # see meaning of status codes here: https://cloufield.github.io/gwaslab/StatusCode/
        missing_from_ref = sumstats.data[GWASLAB_STATUS_COL].str[5:6] == "8"
        print(
            f"Dropping {missing_from_ref.sum()} variants that are missing from the reference"
        )
        sumstats.data = sumstats.data.loc[~missing_from_ref, :]


@frozen
class GWASLabColumnSpecifiers:
    rsid: str | None
    snpid: str | None
    chrom: str
    pos: str
    ea: str
    nea: str
    OR: str | None
    se: str | None
    p: str
    info: str | None
    eaf: str | None = None
    neaf: str | None = None
    beta: str | None = None
    ncase: str | None = None
    ncontrol: str | None = None
    neff: str | None = None
    n: str | None = None
    or_95l: str | None = None
    or_95u: str | None = None


def _get_sumstats(
    x: narwhals.LazyFrame,
    fmt: GwaslabKnownFormat | GWASLabColumnSpecifiers,
    drop_cols: Sequence[str],
) -> gl.Sumstats:
    x = x.drop(drop_cols)
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

    _df_source_task: Task
    _asset_id: AssetId
    basic_check: bool
    genome_build: GenomeBuildMode
    filter_hapmap3: bool = False
    filter_indels: bool = False
    filter_palindromic: bool = False
    exclude_hla: bool = False
    exclude_sexchr: bool = False
    harmonize_options: HarmonizationOptions | None = None
    liftover_to: GenomeBuild | None = None
    fmt: GwaslabKnownFormat | GWASLabColumnSpecifiers = "regenie"
    drop_col_list: Sequence[str] = tuple()
    pre_pipe: DataProcessingPipe = IdentityPipe()

    def __attrs_post_init__(self):
        assert self._source_meta is not None

    @property
    def _source_id(self) -> AssetId:
        return self._df_source_task.meta.asset_id

    @property
    def _source_meta(self) -> FilteredGWASDataMeta | GWASSummaryDataFileMeta:
        meta = self._df_source_task.meta
        assert isinstance(meta, (FilteredGWASDataMeta, GWASSummaryDataFileMeta))
        return meta

    @property
    def meta(self) -> Meta:
        return GWASLabSumStatsMeta(
            short_id=self._asset_id,
            trait=self._source_meta.trait,
            project=self._source_meta.project,
            sub_dir="gwaslab_sumstats",
        )

    @property
    def deps(self) -> list["Task"]:
        return [self._df_source_task]

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> FileAsset:
        df = scan_dataframe_asset(asset=fetch(self._source_id), meta=self._source_meta)
        df = self.pre_pipe.process(df)
        logger.debug("Fetching source dataframe asset...")
        sumstats = _get_sumstats(df, self.fmt, drop_cols=self.drop_col_list)
        transform_spec = GwasLabTransformSpec(
            basic_check=self.basic_check,
            genome_build=self.genome_build,
            filter_hapmap3=self.filter_hapmap3,
            filter_indels=self.filter_indels,
            filter_palindromic=self.filter_palindromic,
            exclude_hla=self.exclude_hla,
            exclude_sexchr=self.exclude_sexchr,
            harmonize_options=self.harmonize_options,
            liftover_to=self.liftover_to,
            fmt=self.fmt,
        )
        sumstats = transform_gwaslab_sumstats(sumstats, spec=transform_spec)
        out_path = scratch_dir / "pickled_sumstats.pickle"
        gl.dump_pickle(sumstats, path=out_path)
        return FileAsset(out_path)


def _sumstats_raise_on_error(sumstats: gl.Sumstats):
    error_status = sumstats.data[GWASLAB_STATUS_COL] == "9999999"
    if error_status.any():
        raise ValueError("GWASLAB Error")
    if len(sumstats.data) == 0:
        raise ValueError("No rows survive GWASLAB quality contro!")


@frozen
class GwasLabTransformSpec:
    basic_check: bool = True
    genome_build: GenomeBuildMode = "infer"
    filter_hapmap3: bool = False
    filter_indels: bool = False
    filter_palindromic: bool = False
    exclude_hla: bool = False
    exclude_sexchr: bool = False
    harmonize_options: HarmonizationOptions | None = None
    liftover_to: GenomeBuild | None = None
    fmt: GwaslabKnownFormat | GWASLabColumnSpecifiers = "regenie"


def transform_gwaslab_sumstats(
    sumstats: gl.Sumstats,
    spec: GwasLabTransformSpec,
    fetch: Fetch | None = None,
) -> gl.Sumstats:
    logger.debug("Running gwas summary statistics through gwaslab pipelines...")
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
    if spec.harmonize_options is not None:
        _do_harmonization(
            sumstats,
            basic_check=(not spec.basic_check),
            options=spec.harmonize_options,
        )

    _sumstats_raise_on_error(sumstats)
    logger.debug(f"Finished gwaslab pipe.  Data has shape {sumstats.data.shape}")
    return sumstats
