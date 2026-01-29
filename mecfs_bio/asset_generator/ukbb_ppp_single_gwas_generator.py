from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.gwas_summary_file_meta import GWASSummaryDataFileMeta
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameTextFormat,
)
from mecfs_bio.build_system.task.concat_frames_in_dir_task import ConcatFramesInDirTask
from mecfs_bio.build_system.task.extract_tar_gzip_task import ExtractTarGzipTask
from mecfs_bio.build_system.task.get_file_from_synapse_task import (
    GetFileFromSynapseTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_region_plot_task_arbitrary_locus import (
    GwasLabRegionPlotTargetLocusTask,
)


@frozen
class UKBBPPPGWASPrep:
    """
    Grouping of tasks for downloading and preprocessing data from a UK Biobank Pharma Proteomics Project GWAS.
    """

    fetch_task: GetFileFromSynapseTask
    untar_task: ExtractTarGzipTask
    stack_task: ConcatFramesInDirTask
    sumstats_37_task: GWASLabCreateSumstatsTask
    plot_task: GwasLabRegionPlotTargetLocusTask


def ubbb_ppp_gwas_prep(
    gene_name: str,
    syn_id: str,
    expected_filename: str,
    region_to_plot_chrom: int,
    region_to_plot_37_pos: int,
):
    """
    Asset generator to generate tasks to download and preprocess data from one of the UK Biobank Pharma Proteomics Project GWAS

    Also plots a region of the downloaded GWAS.
    """
    base_name = gene_name + "_ukbb_ppp"
    trait = "ukbb_ppp"
    fetch_task = GetFileFromSynapseTask(
        meta=GWASSummaryDataFileMeta(
            short_id=AssetId(base_name + "_downloaded_raw_data"),
            trait=trait,
            project=gene_name,
            sub_dir="raw",
            project_path=PurePath(expected_filename),
        ),
        synid=syn_id,
        expected_filename=expected_filename,
    )
    untar_task = ExtractTarGzipTask.create(
        asset_id=base_name + "_untar", source_task=fetch_task, read_mode="r"
    )

    stack_task = ConcatFramesInDirTask.create(
        asset_id=base_name + "_stacked",
        source_dir_task=untar_task,
        path_glob="*.gz",
        read_spec_for_frames=DataFrameReadSpec(DataFrameTextFormat(separator=" ")),
    )
    sumstats_37_task = GWASLabCreateSumstatsTask(
        df_source_task=stack_task,
        asset_id=AssetId(base_name + "_37_sumstats"),
        basic_check=True,
        genome_build="infer",
        liftover_to="19",
        fmt="regenie",
        harmonize_options=None,
    )
    plot_task = GwasLabRegionPlotTargetLocusTask.create(
        asset_id=base_name + "_region_plot",
        sumstats_task=sumstats_37_task,
        vcf_name_for_ld="1kg_eur_hg19",
        chrom=region_to_plot_chrom,
        pos=region_to_plot_37_pos,
    )
    return UKBBPPPGWASPrep(
        fetch_task=fetch_task,
        untar_task=untar_task,
        stack_task=stack_task,
        sumstats_37_task=sumstats_37_task,
        plot_task=plot_task,
    )
