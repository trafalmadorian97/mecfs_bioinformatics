"""
Task generator for fast all-pairs genetic correlation by cross-trait LD score regression.

Each trait is aligned to the shared HapMap3 regression SNPs once (one cached Task per trait), and a
single Task then computes every pair's genetic correlation from the stacked aligned traits, with a
second Task reporting the per-trait heritabilities. Adding a trait therefore prepares only that
trait and reruns the cheap all-pairs step.
"""

from typing import Mapping, Sequence

import narwhals
from attrs import frozen
from frozendict import frozendict

from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.hapmap3_eur_ldsc_variant_index import (
    HAPMAP3_EUR_LDSC_VARIANT_INDEX,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
    ParquetWriteOptions,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.aligned_trait_task import (
    FastCtLdscAlignedTraitTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.alignment_guards import AlignmentGuards
from mecfs_bio.build_system.task.fast_ct_ldsc.all_pairs_rg_task import (
    DEFAULT_RG_BATCH_SIZE,
    FastCtLdscAllPairsRgTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.heritability_task import (
    FastCtLdscHeritabilityTask,
)
from mecfs_bio.build_system.task.fast_ct_ldsc.variant_index_task import (
    FastCtLdscVariantIndexTask,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    SumstatsSource,
)
from mecfs_bio.build_system.task.pipe_dataframe_task import PipeDataFrameTask
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.drop_col_pipe import DropColPipe
from mecfs_bio.build_system.task.pipes.format_numbers_pipe import FormatFloatNumbersPipe
from mecfs_bio.build_system.task.ppp_ldsc.batched_ldsc_h2 import DEFAULT_N_BLOCKS
from mecfs_bio.constants.fast_ct_ldsc_constants import (
    FAST_CT_LDSC_GCOV_COL,
    FAST_CT_LDSC_GCOV_INT_COL,
    FAST_CT_LDSC_GCOV_INT_SE_COL,
    FAST_CT_LDSC_GCOV_SE_COL,
    FAST_CT_LDSC_H2_INT_COL,
    FAST_CT_LDSC_H2_OBS_P1_COL,
    FAST_CT_LDSC_H2_OBS_P2_COL,
    FAST_CT_LDSC_N_SNPS_COL,
    FAST_CT_LDSC_RG_COL,
    FAST_CT_LDSC_RG_P_COL,
    FAST_CT_LDSC_RG_SE_COL,
    FAST_CT_LDSC_RG_Z_COL,
)


@frozen(slots=True)
class FastGeneticCorrTasks:
    aligned_trait_tasks: Mapping[str, FastCtLdscAlignedTraitTask]
    pairs_task: FastCtLdscAllPairsRgTask
    heritability_task: FastCtLdscHeritabilityTask
    pairs_markdown_task: Task
    pairs_display_parquet_task: PipeDataFrameTask

    def terminal_tasks(self) -> list[Task]:
        return [
            self.pairs_task,
            self.heritability_task,
            self.pairs_markdown_task,
            self.pairs_display_parquet_task,
        ]


# Columns of the pairs table left out of the human-facing displays (markdown and display parquet).
_DISPLAY_DROPPED_COLUMNS = [
    FAST_CT_LDSC_GCOV_COL,
    FAST_CT_LDSC_GCOV_SE_COL,
    FAST_CT_LDSC_H2_OBS_P1_COL,
    FAST_CT_LDSC_H2_OBS_P2_COL,
    FAST_CT_LDSC_H2_INT_COL,
]


def fast_genetic_corr_by_ct_ldsc_asset_generator(
    base_name: str,
    sources: Sequence[SumstatsSource],
    variant_index_task: FastCtLdscVariantIndexTask = HAPMAP3_EUR_LDSC_VARIANT_INDEX,
    guards: AlignmentGuards = AlignmentGuards(),
    n_blocks: int = DEFAULT_N_BLOCKS,
    batch_size: int = DEFAULT_RG_BATCH_SIZE,
) -> FastGeneticCorrTasks:
    """Create one aligned-trait Task per source, plus the all-pairs rg, heritability,
    pairs-markdown and pairs display-parquet Tasks. The display parquet holds the same columns as
    the markdown table, unformatted, for the data_table documentation macro. Traits are ordered by lower-cased alias, so p1 < p2 in every pair."""
    ordered = sorted(sources, key=lambda source: source.alias.lower())
    aligned = {
        source.alias: FastCtLdscAlignedTraitTask.create(
            f"{base_name}_aligned_{source.alias.lower()}",
            source=source,
            index_task=variant_index_task,
            guards=guards,
        )
        for source in ordered
    }
    aligned_tasks = list(aligned.values())
    pairs = FastCtLdscAllPairsRgTask.create(
        f"{base_name}_fast_ct_ldsc_pairs",
        index_task=variant_index_task,
        aligned_tasks=aligned_tasks,
        n_blocks=n_blocks,
        batch_size=batch_size,
    )
    heritability = FastCtLdscHeritabilityTask.create(
        f"{base_name}_fast_ct_ldsc_heritability",
        index_task=variant_index_task,
        aligned_tasks=aligned_tasks,
        n_blocks=n_blocks,
    )
    pairs_markdown = ConvertDataFrameToMarkdownTask.create_from_result_table_task(
        asset_id=f"{base_name}_fast_ct_ldsc_pairs_markdown",
        source_task=pairs,
        pipe=CompositePipe(
            [DropColPipe(_DISPLAY_DROPPED_COLUMNS)]
            + [
                FormatFloatNumbersPipe(col=col, format_str=".4g")
                for col in [
                    FAST_CT_LDSC_RG_COL,
                    FAST_CT_LDSC_RG_SE_COL,
                    FAST_CT_LDSC_RG_Z_COL,
                    FAST_CT_LDSC_RG_P_COL,
                    FAST_CT_LDSC_GCOV_INT_COL,
                    FAST_CT_LDSC_GCOV_INT_SE_COL,
                ]
            ]
        ),
    )
    pairs_display_parquet = PipeDataFrameTask.create(
        source_task=pairs,
        asset_id=f"{base_name}_fast_ct_ldsc_pairs_display_parquet",
        pipes=[
            DropColPipe(_DISPLAY_DROPPED_COLUMNS),
            # The data_table macro's parquet reader decodes INT64 as BigInt, which it cannot
            # format.
            CastPipe(
                target_column=FAST_CT_LDSC_N_SNPS_COL,
                type=narwhals.Int32(),
                new_col_name=FAST_CT_LDSC_N_SNPS_COL,
            ),
        ],
        out_format=ParquetOutFormat(
            write_options=ParquetWriteOptions(
                compression="zstd",
                compression_level=22,
                byte_stream_split_floats=True,
            )
        ),
        backend="polars",
    )
    return FastGeneticCorrTasks(
        aligned_trait_tasks=frozendict(aligned),
        pairs_task=pairs,
        heritability_task=heritability,
        pairs_markdown_task=pairs_markdown,
        pairs_display_parquet_task=pairs_display_parquet,
    )
