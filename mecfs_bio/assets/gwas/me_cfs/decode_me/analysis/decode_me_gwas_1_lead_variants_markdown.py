"""
Markdown table of the DecodeME GWAS lead variants, for embedding in the docs in
place of the hand-pasted table.

The gwaslab lead-variants task emits many columns; we select the subset shown in
the write-up.  ``DECODE_ME_GWAS_1_LEAD_VARIANTS`` carries a
``GWASLabLeadVariantsMeta`` (a ``FileMeta`` with a CSV read spec) rather than a
``ResultTableMeta``, so we build the ``MarkdownFileMeta`` directly instead of via
``create_from_result_table_task``.
"""

from pathlib import PurePath

from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_lead_variants import (
    DECODE_ME_GWAS_1_LEAD_VARIANTS,
)
from mecfs_bio.build_system.meta.asset_id import AssetId
from mecfs_bio.build_system.meta.markdown_file_meta import MarkdownFileMeta
from mecfs_bio.build_system.task.convert_dataframe_to_markdown_task import (
    ConvertDataFrameToMarkdownTask,
)
from mecfs_bio.build_system.task.pipes.select_pipe import SelectColPipe

_LEAD_VARIANT_DOC_COLUMNS = [
    "SNPID",
    "CHR",
    "POS",
    "EA",
    "NEA",
    "EAF",
    "BETA",
    "SE",
    "CHISQ",
    "MLOG10P",
    "N",
    "GENE",
]

DECODE_ME_GWAS_1_LEAD_VARIANTS_MARKDOWN = ConvertDataFrameToMarkdownTask(
    meta=MarkdownFileMeta(
        id=AssetId("decode_me_gwas_1_lead_variants_markdown"),
        trait=DECODE_ME_GWAS_1_LEAD_VARIANTS.meta.trait,
        project=DECODE_ME_GWAS_1_LEAD_VARIANTS.meta.project,
        sub_dir=PurePath("analysis/lead_variants"),
    ),
    df_task=DECODE_ME_GWAS_1_LEAD_VARIANTS,
    pipe=SelectColPipe(cols_to_select=_LEAD_VARIANT_DOC_COLUMNS),
)
