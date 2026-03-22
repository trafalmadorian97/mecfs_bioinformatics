"""
This task uses LDSC's wrapper around  the original LDSC to implement cell-type specific analysis
via Stratified LDSC as described in "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types"

See https://github.com/bulik/ldsc/wiki/Cell-type-specific-analyses

and also

https://cloufield.github.io/gwaslab/LDSCinGWASLab/
"""

import shutil
from pathlib import Path, PurePath
from tempfile import TemporaryDirectory

import narwhals
import pandas as pd
import structlog
from attrs import frozen

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwaslab_meta.gwaslab_sumstats_meta import (
    GWASLabSumStatsMeta,
)
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.meta.read_spec.read_sumstats import read_sumstats
from mecfs_bio.build_system.meta.result_table_meta import ResultTableMeta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.wf.base_wf import WF

logger = structlog.get_logger()


@frozen
class CellAnalysisByLDSCTask(Task):
    """
        Task to run cell specific analysis using stratified LDSC
        as described in
    "Heritability enrichment of specifically expressed genes identifies disease-relevant tissues and cell types"

        ref_ld_chr_cts_task: task to generate a directory of cell- or tissue-specific LD scores
        ref_ld_chr_cts_filename:  filename of index into LD scores in directory created by task above
         The official site describes this file as follows:
         "This file has two columns. The first column has a label, for example the name of the cell type in question. The second column has a comma delimited list of LD scores to include when doing the regression for that cell type. For example, in Cahoy.ldcts there are three lines, corresponding to three brain cell types. Each line has two sets of LD scores to include: one is the set of LD scores corresponding to the specifically expressed genes in the cell type, while the second one is a "control" gene set of all genes. The result that will be reported will be the regression coefficient for the first set of LD scores in the list.
         "
        ref_ld_chr_task: task to generate a directory of reference LD scores for the 'baseline' genome annotations.
        ref_ld_chr_inner_pattern: string pattern appended to the path to the above directory to select all chromosome files.  Should be in a special glob format, where the @ character is replaced by chromosome number
        w_ld_chr_task: File containing un-partitioned regression weights for the regression SNPS.  these are typically the hapmap3 variants excluding the HLA/MHC region


    """

    _meta: Meta
    source_sumstats_task: Task
    ref_ld_chr_cts_task: Task
    ref_ld_chr_cts_index_filename: str
    ref_ld_chr_task: Task
    ref_ld_chr_inner_pattern: str
    w_ld_chr_task: Task
    w_ld_chr_inner_pattern: str
    pre_pipe: DataProcessingPipe = IdentityPipe()

    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        return [
            self.source_sumstats_task,
            self.ref_ld_chr_cts_task,
            self.ref_ld_chr_task,
            self.w_ld_chr_task,
        ]

    @property
    def _source_sumstats_id(self) -> AssetId:
        return self.source_sumstats_task.asset_id

    @property
    def _ref_ld_chr_cts_id(self) -> AssetId:
        return self.ref_ld_chr_cts_task.asset_id

    @property
    def _ref_ld_chr_id(self) -> AssetId:
        return self.ref_ld_chr_task.asset_id

    @property
    def _w_ld_chr_id(self) -> AssetId:
        return self.w_ld_chr_task.asset_id

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        sumstats_asset = fetch(self._source_sumstats_id)
        sumstats = read_sumstats(sumstats_asset)
        sumstats.data = (
            self.pre_pipe.process(narwhals.from_native(sumstats.data).lazy())
            .collect()
            .to_pandas()
        )
        ref_ld_chr_cts_asset = fetch(self._ref_ld_chr_cts_id)
        assert isinstance(ref_ld_chr_cts_asset, DirectoryAsset)
        ref_ld_chr_cts_index_path = (
            ref_ld_chr_cts_asset.path / self.ref_ld_chr_cts_index_filename
        )
        assert ref_ld_chr_cts_index_path.is_file()
        with TemporaryDirectory() as tmpdir:
            tmp_path = Path(tmpdir)
            new_loc_chr_cts = tmp_path / "temp_ld_chr_cts"
            shutil.copytree(ref_ld_chr_cts_asset.path, new_loc_chr_cts)
            new_loc_chr_cts_index_path = (
                new_loc_chr_cts / self.ref_ld_chr_cts_index_filename
            )
            assert new_loc_chr_cts_index_path.is_file()
            # there is something of a bug in the index format:
            # the ldsc code assumes that it is
            # being run from within the chr_cts directory, and so treats the paths in the index file
            # as paths relative to the current working directory
            # to run ldsc when the chr_cts directory is not the current directory,
            # we need to prepend a path to all paths in the index file
            prepend_path_to_ldcts_file(
                new_loc_chr_cts_index_path,
                prefix_path=new_loc_chr_cts,
                output_path=new_loc_chr_cts_index_path,
            )

            ref_ld_chr_asset = fetch(self._ref_ld_chr_id)
            assert isinstance(ref_ld_chr_asset, DirectoryAsset)
            ref_ld_chr_full_pattern = (
                ref_ld_chr_asset.path / self.ref_ld_chr_inner_pattern
            )
            # assert ref_ld_chr_full_pattern.is_dir()
            w_ld_chr_asset = fetch(self._w_ld_chr_id)
            assert isinstance(w_ld_chr_asset, DirectoryAsset)
            w_ld_chr_full_pattern = w_ld_chr_asset.path / self.w_ld_chr_inner_pattern
            # assert  w_ld_chr_full_pattern.is_dir()
            sumstats.estimate_h2_cts_by_ldsc(
                ref_ld_chr=str(ref_ld_chr_full_pattern),
                ref_ld_chr_cts=str(new_loc_chr_cts_index_path),
                w_ld_chr=str(w_ld_chr_full_pattern),
            )
            ldsc_h2_cts: pd.DataFrame = sumstats.ldsc_h2_cts
            logger.debug(f"cell type specific s-LDSC results: \n \n {ldsc_h2_cts}\n")

            out_path = scratch_dir / "ldsc_h2_cts.csv"
            ldsc_h2_cts.to_csv(out_path, index=False)
            return FileAsset(out_path)

    @classmethod
    def create(
        cls,
        asset_id: str,
        source_sumstats_task: Task,
        ref_ld_chr_cts_task: Task,
        ref_ld_chr_cts_filename: str,
        ref_ld_chr_task: Task,
        ref_ld_chr_inner_dirname: str,
        w_ld_chr_task: Task,
        w_ld_chr_inner_dirname: str,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
    ):
        sumstats_meta = source_sumstats_task.meta
        assert isinstance(sumstats_meta, GWASLabSumStatsMeta)
        meta = ResultTableMeta(
            id=asset_id,
            trait=sumstats_meta.trait,
            project=sumstats_meta.project,
            sub_dir=PurePath("analysis"),
            extension=".csv",
            read_spec=DataFrameReadSpec(DataFrameTextFormat(separator=",")),
        )
        return cls(
            meta=meta,
            source_sumstats_task=source_sumstats_task,
            ref_ld_chr_cts_task=ref_ld_chr_cts_task,
            ref_ld_chr_cts_index_filename=ref_ld_chr_cts_filename,
            ref_ld_chr_task=ref_ld_chr_task,
            ref_ld_chr_inner_pattern=ref_ld_chr_inner_dirname,
            w_ld_chr_task=w_ld_chr_task,
            w_ld_chr_inner_pattern=w_ld_chr_inner_dirname,
            pre_pipe=pre_pipe,
        )


@frozen
class LDCTSFileEntry:
    label: str
    target_gene_data_path: Path
    control_gene_data_path: Path

    def __attrs_post_init__(self):
        forbidden_chars = ["\t", " ", ","]
        for char in forbidden_chars:
            assert char not in self.label
            assert char not in str(self.target_gene_data_path)
            assert char not in str(self.control_gene_data_path)

    def to_text(self) -> str:
        return (
            f"{self.label}\t{self.target_gene_data_path},{self.control_gene_data_path}"
        )

    def prepend_to_paths(self, prefix: Path) -> "LDCTSFileEntry":
        return LDCTSFileEntry(
            label=self.label,
            target_gene_data_path=prefix / self.target_gene_data_path,
            control_gene_data_path=prefix / self.control_gene_data_path,
        )

    @classmethod
    def from_text(cls, text: str) -> "LDCTSFileEntry":
        label, rest = text.split("\t", maxsplit=1)
        target, control = rest.split(",", maxsplit=1)
        return cls(
            label=label,
            target_gene_data_path=Path(target),
            control_gene_data_path=Path(control),
        )


def ldcts_raw_to_processed(contents: str) -> list[LDCTSFileEntry]:
    return [LDCTSFileEntry.from_text(line) for line in contents.splitlines()]


def ldcts_processed_to_raw(processed: list[LDCTSFileEntry]) -> str:
    return "\n".join([item.to_text() for item in processed])


def prepend_path_to_ldcts_file(input_path: Path, prefix_path: Path, output_path: Path):
    with open(input_path) as f:
        contents = f.read()
    structured_contents = ldcts_raw_to_processed(contents)
    prefixed = [item.prepend_to_paths(prefix_path) for item in structured_contents]
    with open(output_path, "w") as f:
        f.write(ldcts_processed_to_raw(prefixed))
