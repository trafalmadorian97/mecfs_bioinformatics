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
from mecfs_bio.build_system.task.gwaslab.gwaslab_util import (
    gwaslab_download_ref_if_missing,
)
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

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
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_STATUS_COL,
    GwaslabKnownFormat,
)

GenomeBuildMode = Literal["infer", "19", "38"]
GenomeBuild = Literal["19", "38"]


@frozen
class GWASLabVCFRef:
    name: str
    ref_alt_freq: str
    extra_downloads: Sequence[str] = tuple()


@frozen
class HarmonizationOptions:
    """
    Options for the call to GWASLab's harmonize function.


    gwaslab's harmonization function changes the status codes in the STATUS column.
    These status codes are described here: https://cloufield.github.io/gwaslab/StatusCode/

    Below I explain some points that were initially not clear to me from the gwaslab documentation.

    Two reference files are used in harmonization:

    - A VCF file (ref_infer).  This is basically a table of genetic variants.

       In some cases, this table is in dbSNP VCF format.  In this case, each row describes a given genetic variant.
       Sometimes this description includes allele frequency.

       In other cases, (such as when using a thousand genomes reference data) this table is in genotype VCF format.
       In this case, the rows of the VCF file correspond to variants, and the columns correspond to individuals (from the
       thousand genomes project, for example).  For each individual and each variant, the table tells us whether that individual has that variant.
       Variant frequency information can be calculated from this individual-level genome data.


    - A FASTA file (ref_seq).  This is a consensus human genome sequence.  Here is an example of some rows from the hg19 FASTA file:

        TAAGTTTTGTCTGGTAATAAAGGTATATTTTCAAAAGAGAGGTAAATAGA
        TCCACATACTGTGGAGGGAATAAAATACTTTTTGAAAAACAAACAACAAG
        TTGGATTTTTAGACACATAGAAATTGAATATGTACATTTATAAATATTTT
        TGGATTGAACTATTTCAAAATTATACCATAAAATAACTTGTAAAAATGTA
        GGCAAAATGTATATAATTATGGCATGAGGTATGCAACTTTAGGCAAGGAA
        GCAAAAGCAGAAACCATGAAAAAAGTCTAAATTTTACCATATTGAATTTA
        AATTTTCAAAAACAAAAATAAAGACAAAGTGGGAAAAATATGTATGCTTC
        ATGTGTGACAAGCCACTGATACCTATTAAATATGAAGAATATTATAAATC
        ATATCAATAACCACAACATTCAAGCTGTCAGTTTGAATAGACaatgtaaa
        tgacaaaactacatactcaacaagataacagcaaaccagcttcgacagca
        cgttaaaggggtcatacaacataatcgagtagaatttatctctgagatgc
        aagaatggttcaaaatatggaaaccaataaatgtgatatgccacactaac
        agaataaaaaataaaaatcatattatcatctcaatagatgcagaaaaagc
        attaacaaaagtaaacattctttcataataagacatcagataaaacaaat
        taggaatagaaggaatgtaccgcaacacaataaaggccatatataacaag
        cccacagctaacatcataatagtaaaatcatcacactggtaaaaaaaatg





    gwaslab uses these two reference files to harmonize summary statistics.
    These two reference files each affect a different digit of the gwaslab STATUS code column.

    - Digit 7 of the status code is determined by the ability of gwaslab to find the variant in the reference VCF (ref_infer) :
      see here: https://github.com/Cloufield/gwaslab/blob/d639b67c5264b1ac7ec89e284e638f2c8454ac48/src/gwaslab/hm/hm_harmonize_sumstats.py#L1521-L1530
      Values of 7 or 8 here mean that the variant is palindromic, and the database could not be used to disambiguate the strand of the variant, or the variant was not found in the database
    - Digit 6 of the status code is instead determined by the ability of the gwaslab to find the variant in the reference genome build FASTA file (ref_seq)
      see here: https://github.com/Cloufield/gwaslab/blob/d639b67c5264b1ac7ec89e284e638f2c8454ac48/src/gwaslab/hm/hm_harmonize_sumstats.py#L968-L975
      a value of 8 means a failure to find the variant in the FASTA file.

    Set drop_missing_from_ref_seq to drop based on digit 6.
    Set drop_missing_from_ref_infer to drop based on digit 7.

    """

    ref_infer: GWASLabVCFRef
    ref_seq: str
    cores: int
    check_ref_files: bool
    drop_missing_from_ref_seq: bool
    drop_missing_from_ref_infer_or_ambiguous: bool = True


def _do_harmonization(
    sumstats: gl.Sumstats, basic_check: bool, options: HarmonizationOptions
):
    if options.check_ref_files:
        gwaslab_download_ref_if_missing(options.ref_infer.name)
        gwaslab_download_ref_if_missing(options.ref_seq)
        for extra in options.ref_infer.extra_downloads:
            gwaslab.download_ref(name=extra, overwrite=False)

    sumstats.harmonize(
        basic_check=basic_check,
        n_cores=options.cores,
        ref_seq=gl.get_path(options.ref_seq),
        ref_infer=gl.get_path(options.ref_infer.name),
        ref_alt_freq=options.ref_infer.ref_alt_freq,
    )
    if options.drop_missing_from_ref_seq:
        # see meaning of status codes here: https://cloufield.github.io/gwaslab/StatusCode/
        missing_from_ref_seq = sumstats.data[GWASLAB_STATUS_COL].str[5:6] == "8"
        logger.debug(
            f"Dropping {missing_from_ref_seq.sum()} variants that are missing from the sequence FASTA reference"
        )
        sumstats.data = sumstats.data.loc[~missing_from_ref_seq, :]
    if options.drop_missing_from_ref_infer_or_ambiguous:
        missing_from_ref_infer = (sumstats.data[GWASLAB_STATUS_COL].str[6] == "8") | (
            sumstats.data[GWASLAB_STATUS_COL].str[6] == "7"
        )
        logger.debug(
            f"Dropping {missing_from_ref_infer.sum()} variants that are missing from the VCF reference"
        )
        sumstats.data = sumstats.data.loc[~missing_from_ref_infer, :]


@frozen
class GWASLabColumnSpecifiers:
    rsid: str | None
    snpid: str | None
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

    def get_selection_pipe(self) -> SelectColPipe:
        fields = attrs.asdict(self)
        cols = []
        for v in fields.values():
            if v is not None:
                cols.append(v)
        return SelectColPipe(cols)


ValidGwaslabFormat = GwaslabKnownFormat | GWASLabColumnSpecifiers


def _get_sumstats(
    x: narwhals.LazyFrame,
    fmt: ValidGwaslabFormat,
    drop_cols: Sequence[str],
) -> gl.Sumstats:
    if isinstance(fmt, GWASLabColumnSpecifiers):
        x = fmt.get_selection_pipe().process(x)
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
            chisq=fmt.chi_sq,
            mlog10p=fmt.mlog10p,
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
    def _source_meta(
        self,
    ) -> FilteredGWASDataMeta | GWASSummaryDataFileMeta | ReferenceFileMeta:
        meta = self._df_source_task.meta
        assert isinstance(
            meta, (FilteredGWASDataMeta, GWASSummaryDataFileMeta, ReferenceFileMeta)
        )
        return meta

    @property
    def meta(self) -> Meta:
        if isinstance(self._source_meta, ReferenceFileMeta):
            return GWASLabSumStatsMeta(
                id=self._asset_id,
                trait="reference_data_gwas",
                project=self._source_meta.group,
            )
        return GWASLabSumStatsMeta(
            id=self._asset_id,
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
    harmonize_options: HarmonizationOptions | None = None
    liftover_to: GenomeBuild | None = None


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
