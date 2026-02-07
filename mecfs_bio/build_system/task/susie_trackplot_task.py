import narwhals
import structlog
import polars as pl
import pandas as pd
import tempfile
from pathlib import Path
import configparser
from attrs import frozen
from rich.pretty import pprint
from statsmodels.tsa.tests.test_bds import data_file

from mecfs_bio.build_system.asset.base_asset import Asset
from mecfs_bio.build_system.asset.directory_asset import DirectoryAsset
from mecfs_bio.build_system.asset.file_asset import FileAsset
from mecfs_bio.build_system.meta.meta import Meta
from mecfs_bio.build_system.rebuilder.fetch.base_fetch import Fetch
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_mlog10p_pipe import ComputeMlog10pIfNeededPipe
from mecfs_bio.build_system.task.pipes.compute_p_from_beta_se import ComputePFromBetaSEPipeIfNeeded
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import COMBINED_CS_FILENAME, FILTERED_GWAS_FILENAME, \
    CS_COLUMN, PIP_COLUMN
from mecfs_bio.build_system.wf.base_wf import WF
from mecfs_bio.constants.gwaslab_constants import GWASLAB_CHROM_COL, GWASLAB_POS_COL, GWASLAB_MLOG10P_COL
from mecfs_bio.util.subproc.run_command import execute_command
logger = structlog.get_logger()

# Try: All credible sets in one track.

@frozen
class  EnsemblGeneInfoSource:
    source_task: Task

GeneInfoSource = EnsemblGeneInfoSource

PLOT_CONFIG_FILENAME = "plot_config.txt"

PLOT_FILENAME = "track_plot.png"


GWAS_BEDGRAPH_FILENAME = "gwas.bedgraph"

@frozen
class RegionSelectOverride:
    chrom: int
    start: int
    end: int

@frozen
class RegionSelectDefault:
    pass
RegionSelect = RegionSelectOverride|RegionSelectDefault




@frozen
class SusieTrackPlotTask(Task):
    _meta: Meta
    susie_task: Task
    gene_info_source: EnsemblGeneInfoSource|None
    region_mode: RegionSelect


    @property
    def meta(self) -> Meta:
        return self._meta

    @property
    def deps(self) -> list["Task"]:
        deps =  [self.susie_task]
        if self.gene_info_source is not None:
            deps.append(self.gene_info_source.source_task)
        return deps

    def execute(self, scratch_dir: Path, fetch: Fetch, wf: WF) -> Asset:
        susie_asset = fetch(self.susie_task.asset_id)
        assert isinstance(susie_asset, DirectoryAsset)
        susie_output_path = susie_asset.path


        chrom, start, end = get_region(self.region_mode, susie_output_path=susie_output_path)

        plot_config = {}

        plot_config = plot_config | get_gwas_tracks(
            susie_output_path=susie_output_path,
            scratch_dir=scratch_dir,
        )

        plot_config = plot_config | get_susie_tracks(susie_dir_path=susie_output_path, scratch_dir=scratch_dir)

        gene_track_file = prepare_gene_info_source(src=self.gene_info_source, fetch=fetch, target_dir=scratch_dir)
        plot_config = plot_config | get_gene_track(gene_track_file)

        plot_config = plot_config |get_x_axis(

        )



        config_path = scratch_dir / PLOT_CONFIG_FILENAME

        pprint(plot_config)
        parser = config_to_conifgparser(plot_config)
        with open(config_path,"w") as f:
            parser.write(f)

        plot_path = scratch_dir/PLOT_FILENAME

        execute_command([
           "pyGenomeTracks",
            "--tracks",
            str(config_path),
            "--region",
            f"{CHROM_PREFIX}{chrom}:{start}-{end}",
            "-o",
            str(plot_path)
            ]
        )
        logger.debug(f"Wrote plot to tracks plot {plot_path}")

        import pdb; pdb.set_trace()

        return DirectoryAsset(scratch_dir)


def get_gwas_tracks(
        susie_output_path: Path,
        scratch_dir: Path,
        max_mlog10p:int=200
        )-> dict:
    filtered_gwas = pl.read_parquet(susie_output_path/FILTERED_GWAS_FILENAME)
    pipe = CompositePipe(
        [

            ComputePFromBetaSEPipeIfNeeded(),
            ComputeMlog10pIfNeededPipe()
        ]
    )
    filtered_gwas = pipe.process(narwhals.from_native(filtered_gwas).lazy()).collect().to_polars()
    filtered_gwas = add_labeled_chrom_col(filtered_gwas)
    filtered_gwas =add_end_col(filtered_gwas)
    filtered_gwas = filtered_gwas.with_columns(
        pl.min_horizontal(pl.lit(max_mlog10p), pl.col(GWASLAB_MLOG10P_COL)).alias(GWASLAB_MLOG10P_COL),
    )
    out_path = scratch_dir/GWAS_BEDGRAPH_FILENAME
    filtered_gwas.select([LABELED_CHROM_COL, GWASLAB_POS_COL, END_COL, GWASLAB_MLOG10P_COL]).write_csv(out_path, separator = "\t", include_header = False
                                                                                                       )
    result = {
        "gwas_track":{
            "file": str(out_path),
            "title": f"GWAS -log_10(p)",
            "color": "green",
            "height": 4,
            "file_type": "bedgraph",
            "min_value": 2,
            # "max_value": 1,
            # "y_axis_values": "0, 1"
        },
        "gwas_spacer":{
            "file_type": "spacer"
        }
    }
    return result


def get_x_axis()->dict:
    return {
        "x-axis":{
        }
    }

def config_to_conifgparser(config: dict) -> configparser.ConfigParser:
    parser = configparser.ConfigParser()
    for key, value in config.items():
        parser[key] = value
    return parser


def get_region(mode:RegionSelect, susie_output_path: Path)->tuple[int, int, int]:
    if isinstance(mode, RegionSelectOverride):
        return mode.chrom, mode.start, mode.end
    df= pl.read_parquet(susie_output_path/FILTERED_GWAS_FILENAME)
    assert df[GWASLAB_CHROM_COL].n_unique()==1
    return df[GWASLAB_CHROM_COL][0], max(df[GWASLAB_POS_COL].min()-1,0), (df[GWASLAB_POS_COL].max() + 2)



def prepare_gene_info_source(src: GeneInfoSource|None, fetch:Fetch, target_dir: Path)->Path|None:
    if src is None:
        return None
    assert isinstance(src, EnsemblGeneInfoSource)
    src_asset = fetch(src.source_task.asset_id)
    assert isinstance(src_asset, FileAsset)

    df = pd.read_csv(src_asset.path,sep="\s+",
                     header=None, names=["ensembl", GWASLAB_CHROM_COL, "start", "end","strand", "name"])
    df = add_labeled_chrom_col(pl.from_pandas(df)).to_pandas()
    out_path = target_dir/"gene_track_file.txt"
    df["intensity"]=100
    df.loc[:,[LABELED_CHROM_COL,"start","end","name", "intensity","strand"]].to_csv(out_path, index=False, sep="\t", header=False)
    return out_path

def get_gene_track(
       track_file_path: Path|None
    ) -> dict:
    if track_file_path is None:
        return {}
    return {
        "genes":{
            "file": str(track_file_path),
            "title":"Genes",
            "height" : 3,
    "fontsize" : 10,
    "file_type" : "bed",
    "gene_rows" : 10
        },
        "genes_spacer":{
            "file_type":"spacer"
        }
    }

def get_susie_tracks(susie_dir_path: Path, scratch_dir: Path)->dict:
    result = {}
    susie_cs_data = pl.read_parquet(susie_dir_path/COMBINED_CS_FILENAME)
    for i,(name_tuple, data) in enumerate(susie_cs_data.group_by(CS_COLUMN)):
        name = name_tuple[0]
        data =   data.with_columns(  (pl.col(GWASLAB_POS_COL)+1).alias("end")  )
        data = add_labeled_chrom_col(data)
        out_path = scratch_dir/ f"{name}_data.bedmap"
        data.select([LABELED_CHROM_COL, GWASLAB_POS_COL, "end", PIP_COLUMN]).write_csv(out_path, separator="\t", include_header=False)
        result[name]={
            "file": str(out_path),
            "title":f"CS {i+1} PIP",
            "color":"red",
            "height":3,
            "file_type":"bedgraph",
            "min_value":0,
            "max_value":1,
            "y_axis_values" : "0, 1"

        # "type":"points:4"
        }
        result[name+"_space"]={
            "file_type" : "spacer"
        }
    return result

CHROM_PREFIX = "chr"
LABELED_CHROM_COL = "labeled_chrom"

END_COL ="end"

def add_labeled_chrom_col(df: pl.DataFrame) -> pl.DataFrame:
    return df.with_columns(
        (pl.lit(CHROM_PREFIX)+ pl.col(GWASLAB_CHROM_COL).cast(pl.String()) ).alias(LABELED_CHROM_COL),
    )

def add_end_col(df:pl.DataFrame) -> pl.DataFrame:
    return df.with_columns( (pl.col(GWASLAB_POS_COL)+1).alias(END_COL), )
