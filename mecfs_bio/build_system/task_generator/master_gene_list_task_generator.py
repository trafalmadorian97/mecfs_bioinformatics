from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.task.annotate_genes_from_reference_task import (
    AnnotateGenesFromReferenceTask,
    GeneInfoReferenceTasks,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.combine_gene_lists_task import (
    ENSEMBL_ID_LABEL,
    CombineGeneListsTask,
    SrcGeneList,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe
from mecfs_bio.constants.gene_info_constants import (
    GENE_INFO_NAME_COL,
    GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
    GENE_INFO_SYMBOL_COL,
    GENE_INFO_UNIPROT_DESCRIPTION_COL,
    GENE_INFO_UNIPROT_ID_COL,
)


@frozen(slots=True)
class MasterGeneListTasks:
    combine_task: Task
    gene_info_task: Task
    markdown_task: Task

    def terminal_tasks(self) -> list[Task]:
        return [self.markdown_task]


def generate_master_gene_list_tasks(
    base_name: str,
    gene_sources: Sequence[SrcGeneList],
    gene_info_reference: GeneInfoReferenceTasks,
) -> MasterGeneListTasks:
    combine_task = CombineGeneListsTask.create(
        asset_id=base_name + "_combine_gene_lists",
        src_gene_lists=gene_sources,
        out_format=ParquetOutFormat(),
    )
    gene_info_task = AnnotateGenesFromReferenceTask.create(
        asset_id=base_name + "_master_gene_list_with_gene_info",
        source_df_task=combine_task,
        ensembl_id_col=ENSEMBL_ID_LABEL,
        reference=gene_info_reference,
        post_pipe=SortPipe([ENSEMBL_ID_LABEL]),
        out_format=ParquetOutFormat(),
    )
    markdown_task = ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        gene_info_task,
        base_name + "_combined_gene_list_markdown",
        pipe=CompositePipe(
            [
                SelectColPipe(
                    [
                        ENSEMBL_ID_LABEL,
                        "sources",
                        GENE_INFO_SYMBOL_COL,
                        GENE_INFO_NAME_COL,
                        GENE_INFO_UNIPROT_ID_COL,
                        GENE_INFO_SUBCELLULAR_LOCALISATION_COL,
                        GENE_INFO_UNIPROT_DESCRIPTION_COL,
                    ]
                ),
                DropNullsPipe([ENSEMBL_ID_LABEL]),
            ]
        ),
    )
    return MasterGeneListTasks(
        combine_task=combine_task,
        gene_info_task=gene_info_task,
        markdown_task=markdown_task,
    )
