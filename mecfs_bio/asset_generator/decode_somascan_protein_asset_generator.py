"""
Asset generator for the deCODE SomaScan database: one BuildDecodeSomascanProteinTask per assay of
a group's manifest, aligned onto that group's shared variant index.

Asset ids are namespaced by group, index, gene symbol and seq id; the gene symbol alone is not
unique across assays, and the seq id's "-" is replaced because it collides with the asset id
delimiter.
"""

from attrs import frozen

from mecfs_bio.asset_generator.decode_somascan_manifest import (
    load_decode_somascan_proteins,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.decode_database.build_decode_somascan_protein_task import (
    DECODE_SLIM_MEMBER,
    BuildDecodeSomascanProteinTask,
)
from mecfs_bio.build_system.task.protein_coloc.protein_hit_table_task import (
    HitTableProtein,
)
from mecfs_bio.constants.decode_somascan_constants import DecodeSomascanGroup
from mecfs_bio.constants.protein_coloc_constants import ProteinId


@frozen(slots=True)
class DecodeSomascanTaskCollection:
    group: DecodeSomascanGroup
    protein_tasks: tuple[BuildDecodeSomascanProteinTask, ...]

    def terminal_tasks(self) -> list[Task]:
        return list(self.protein_tasks)

    def hit_table_proteins(self) -> tuple[HitTableProtein, ...]:
        """Each protein's slim file, labelled by seq id and gene symbol."""
        return tuple(
            HitTableProtein(
                task=task,
                protein_id=ProteinId(task.protein.seq_id),
                gene=task.protein.gene_symbol,
                slim_member=DECODE_SLIM_MEMBER,
            )
            for task in self.protein_tasks
        )


def decode_protein_asset_id(
    group: DecodeSomascanGroup, index_name: str, gene_symbol: str, seq_id: str
) -> str:
    return (
        f"decode_somascan_{group}_{index_name}_{gene_symbol}_{seq_id.replace('-', '_')}"
    )


def generate_decode_somascan_protein_tasks(
    group: DecodeSomascanGroup, index_task: Task, index_name: str
) -> DecodeSomascanTaskCollection:
    """
    index_task: the group's ConstructDecodeVariantIndexTask.
    index_name: a short label for the index (e.g. hapmap_3), namespacing asset ids and paths.
    """
    return DecodeSomascanTaskCollection(
        group=group,
        protein_tasks=tuple(
            BuildDecodeSomascanProteinTask.create(
                index_task=index_task,
                protein=protein,
                asset_id=decode_protein_asset_id(
                    group, index_name, protein.gene_symbol, protein.seq_id
                ),
                index_name=index_name,
                group=group,
            )
            for protein in load_decode_somascan_proteins(group)
        ),
    )
