from typing import Sequence

from attrs import frozen

from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.combine_gene_lists_task import (
    ENSEMBL_ID_LABEL,
    CombineGeneListsTask,
    SrcGeneList,
)
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.fetch_gget_info_task import FetchGGetInfoTask
from mecfs_bio.build_system.task.pipe_dataframe_task import ParquetOutFormat
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_null_pipe import DropNullsPipe
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe
from mecfs_bio.build_system.task.pipes.sort_pipe import SortPipe


@frozen
class MasterGeneListTasks:
    combine_task: Task
    gget_label_task: Task
    markdown_task: Task

    def terminal_tasks(self) -> list[Task]:
        return [self.markdown_task]


def generate_master_gene_list_tasks(
    base_name: str, gene_sources: Sequence[SrcGeneList], max_genes_to_use: int = 100
) -> MasterGeneListTasks:
    combine_task = CombineGeneListsTask.create(
        asset_id=base_name + "_combine_gene_lists",
        src_gene_lists=gene_sources,
        out_format=ParquetOutFormat(),
    )
    gget_label_task = FetchGGetInfoTask.create(
        asset_id=base_name + "_master_gene_list_with_gget",
        source_df_task=combine_task,
        ensembl_id_col=ENSEMBL_ID_LABEL,
        post_pipe=SortPipe([ENSEMBL_ID_LABEL]),
        out_format=ParquetOutFormat(),
        genes_to_use=max_genes_to_use,
    )
    markdown_task = ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        gget_label_task,
        base_name + "_combined_gene_list_markdown",
        pipe=CompositePipe(
            [
                SelectColPipe(
                    [
                        ENSEMBL_ID_LABEL,
                        "sources",
                        "uniprot_id",
                        "primary_gene_name",
                        "subcellular_localisation",
                        "ensembl_description",
                        "uniprot_description",
                        "ncbi_description",
                    ]
                ),
                DropNullsPipe([ENSEMBL_ID_LABEL]),
            ]
        ),
    )
    return MasterGeneListTasks(
        combine_task=combine_task,
        gget_label_task=gget_label_task,
        markdown_task=markdown_task,
    )
