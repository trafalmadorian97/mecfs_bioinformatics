from pathlib import PurePath

from attrs import frozen

from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.reference.schemas.chrom_rename_rules import (
    CHROM_RENAME_RULES,
)
from mecfs_bio.build_system.reference.schemas.hg19_snp151_schema_valid_choms import (
    HG19_SNP151_VALID_CHROMS,
)
from mecfs_bio.build_system.task.assign_rsids_via_snp151_task import (
    AssignRSIDSToSNPsViaSNP151Task,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.fdr_multiple_testing_table_task import (
    MultipleTestingTableTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_constants import GwaslabKnownFormat
from mecfs_bio.build_system.task.gwaslab.gwaslab_create_sumstats_task import (
    GenomeBuildMode,
    GWASLabColumnSpecifiers,
    GWASLabCreateSumstatsTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.magma.magma_annotate_task import MagmaAnnotateTask
from mecfs_bio.build_system.task.magma.magma_gene_analysis_task import (
    MagmaGeneAnalysisTask,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    MagmaGeneSetAnalysisTask,
    ModelParams,
)
from mecfs_bio.build_system.task.magma.magma_plot_gene_set_result import (
    MAGMAPlotGeneSetResult,
)
from mecfs_bio.build_system.task.magma.magma_snp_location_task import MagmaSNPFileTask
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe


@frozen
class StandardMagmaTaskGenerator:
    """
    A Task Generator that generates all the Tasks required to perform the main steps of the MAGMA workflow, up to plotting the results of tissue-specific gene-set analysis
    """

    p_value_task: Task
    snp_loc_task: Task
    annotations_task: Task
    gene_analysis_task: Task
    gene_set_analysis_task: Task
    bar_plot_task: Task
    filtered_gene_analysis_task: Task
    labeled_filtered_gene_analysis_task: Task | None = None

    @classmethod
    def create(
        cls,
        gwas_parquet_with_rsids_task: Task,
        magma_binary_task: Task,
        gene_loc_file_task: Task,
        magma_ld_ref_task: Task,
        tissue_expression_gene_set_task: Task,
        base_name: str,
        sample_size: int,
        ld_ref_file_stem: str = "g1000_eur",
        gene_thesaurus_task: Task | None = None,
    ):
        p_value_task = (
            MagmaSNPFileTask.create_for_magma_snp_p_value_file_compute_if_needed(
                gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
                asset_id=base_name + "_magma_snp_p_values",
            )
        )
        snp_loc_task = MagmaSNPFileTask.create_for_magma_snp_pos_file(
            gwas_parquet_with_rsids_task=gwas_parquet_with_rsids_task,
            asset_id=base_name + "_magma_snp_locs",
        )
        annotations_task = MagmaAnnotateTask.create(
            asset_id=base_name + "_magma_annotations",
            magma_binary_task=magma_binary_task,
            snp_loc_file_task=snp_loc_task,
            gene_loc_file_task=gene_loc_file_task,
        )
        gene_analysis_task = MagmaGeneAnalysisTask.create(
            asset_id=base_name + "_magma_gene_analysis",
            magma_annotation_task=annotations_task,
            magma_p_value_task=p_value_task,
            magma_binary_task=magma_binary_task,
            magma_ld_ref_task=magma_ld_ref_task,
            ld_ref_file_stem=ld_ref_file_stem,
            sample_size=sample_size,
        )
        tissue_gene_set_analysis = MagmaGeneSetAnalysisTask.create(
            asset_id=base_name + "_magma_tissue_gene_set_analysis",
            magma_gene_analysis_task=gene_analysis_task,
            magma_binary_task=magma_binary_task,
            gene_set_task=tissue_expression_gene_set_task,
            set_or_covar="covar",
            sub_dir=PurePath("analysis"),
            model_params=ModelParams(
                direction_covar="greater", condition_hide=["Average"]
            ),
        )
        bar_plot_task = MAGMAPlotGeneSetResult.create(
            gene_set_analysis_task=tissue_gene_set_analysis,
            asset_id=base_name + "_magma_bar_plot",
        )

        filtered_gene_task = (
            MultipleTestingTableTask.create_from_magma_gene_analysis_task(
                asset_id=base_name + "_magma_gene_analysis_filtered",
                alpha=0.01,
                method="fdr_bh",
                source_task=gene_analysis_task,
            )
        )
        if gene_thesaurus_task is not None:
            labeled_filtered_gene_task = JoinDataFramesTask.create_from_result_df(
                asset_id=base_name + "_magma_gene_analysis_filtered_labeled",
                result_df_task=filtered_gene_task,
                reference_df_task=gene_thesaurus_task,
                how="left",
                left_on=["GENE"],
                right_on=["Gene stable ID"],
            )
        else:
            labeled_filtered_gene_task = None

        return cls(
            p_value_task=p_value_task,
            snp_loc_task=snp_loc_task,
            annotations_task=annotations_task,
            gene_analysis_task=gene_analysis_task,
            gene_set_analysis_task=tissue_gene_set_analysis,
            bar_plot_task=bar_plot_task,
            filtered_gene_analysis_task=filtered_gene_task,
            labeled_filtered_gene_analysis_task=labeled_filtered_gene_task,
        )


@frozen
class MagmaTaskGeneratorFromRaw:
    """
    As above, but starts from raw data.
    Assumes raw data already includes RSIDs
    """

    sumstats_task: Task
    parquet_file_task: Task
    inner: StandardMagmaTaskGenerator

    @classmethod
    def create(
        cls,
        raw_gwas_data_task: Task,
        fmt: GwaslabKnownFormat | GWASLabColumnSpecifiers,
        magma_binary_task: Task,
        gene_loc_file_task: Task,
        magma_ld_ref_task: Task,
        tissue_expression_gene_set_task: Task,
        base_name: str,
        sample_size: int,
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        post_pipe: DataProcessingPipe = IdentityPipe(),
        ld_ref_file_stem: str = "g1000_eur",
        genome_build: GenomeBuildMode = "infer",
        gene_thesaurus_task: Task | None = None,
    ):
        sumstats_task = GWASLabCreateSumstatsTask(
            df_source_task=raw_gwas_data_task,
            asset_id=AssetId(base_name + "_sumstats_37"),
            basic_check=True,
            genome_build=genome_build,
            liftover_to="19",
            fmt=fmt,
            pre_pipe=pre_pipe,
        )
        parquet_file_task = GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=sumstats_task,
            asset_id=base_name + "_parquet_table_from_sumstats",
            sub_dir="processed",
            pipe=post_pipe,
        )
        return cls(
            sumstats_task=sumstats_task,
            parquet_file_task=parquet_file_task,
            inner=StandardMagmaTaskGenerator.create(
                gwas_parquet_with_rsids_task=parquet_file_task,
                magma_binary_task=magma_binary_task,
                gene_loc_file_task=gene_loc_file_task,
                magma_ld_ref_task=magma_ld_ref_task,
                tissue_expression_gene_set_task=tissue_expression_gene_set_task,
                base_name=base_name,
                sample_size=sample_size,
                ld_ref_file_stem=ld_ref_file_stem,
                gene_thesaurus_task=gene_thesaurus_task,
            ),
        )


@frozen
class MagmaTaskGeneratorFromRawCompute37RSIDs:
    """
    As above, but assumes RSIDs are either not present,
    or present but invalid because the data needs to be lifter over to assembly 37
    """

    sumstats_task: Task
    parquet_file_task: Task
    assign_rsids_task: Task
    inner: StandardMagmaTaskGenerator

    @classmethod
    def create(
        cls,
        raw_gwas_data_task: Task,
        fmt: GwaslabKnownFormat | GWASLabColumnSpecifiers,
        magma_binary_task: Task,
        gene_loc_file_task: Task,
        magma_ld_ref_task: Task,
        snp151_database_file_task: Task,
        tissue_expression_gene_set_task: Task,
        base_name: str,
        sample_size: int,
        pipe: DataProcessingPipe = IdentityPipe(),
        pre_pipe: DataProcessingPipe = IdentityPipe(),
        ld_ref_file_stem: str = "g1000_eur",
        genome_build: GenomeBuildMode = "infer",
    ):
        sumstats_task = GWASLabCreateSumstatsTask(
            df_source_task=raw_gwas_data_task,
            asset_id=AssetId(base_name + "_sumstats_37"),
            basic_check=True,
            genome_build=genome_build,
            liftover_to="19",
            fmt=fmt,
            pre_pipe=pre_pipe,
        )
        parquet_file_task = GwasLabSumstatsToTableTask.create_from_source_task(
            source_tsk=sumstats_task,
            asset_id=base_name + "_parquet_table_from_sumstats",
            sub_dir="processed",
            pipe=pipe,
        )
        assign_rsids_37_task = AssignRSIDSToSNPsViaSNP151Task.create(
            snp151_database_file_task=snp151_database_file_task,
            raw_snp_data_task=parquet_file_task,
            asset_id=base_name + "_assign_rsids",
            valid_chroms=HG19_SNP151_VALID_CHROMS,
            chrom_replace_rules=CHROM_RENAME_RULES,
        )

        return cls(
            sumstats_task=sumstats_task,
            parquet_file_task=parquet_file_task,
            inner=StandardMagmaTaskGenerator.create(
                gwas_parquet_with_rsids_task=assign_rsids_37_task,
                magma_binary_task=magma_binary_task,
                gene_loc_file_task=gene_loc_file_task,
                magma_ld_ref_task=magma_ld_ref_task,
                tissue_expression_gene_set_task=tissue_expression_gene_set_task,
                base_name=base_name,
                sample_size=sample_size,
                ld_ref_file_stem=ld_ref_file_stem,
            ),
            assign_rsids_task=assign_rsids_37_task,
        )
