from pathlib import PurePath

from mecfs_bio.assets.gwas.schizophrenia.pgc2022.analysis.magma.pgc2022_sch_magma_human_brain_atlas import (
    MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR,
)
from mecfs_bio.build_system.meta.read_spec.dataframe_read_spec import (
    DataFrameReadSpec,
    DataFrameWhiteSpaceSepTextFormat,
)
from mecfs_bio.build_system.task.copy_file_from_directory_task import (
    CopyFileFromDirectoryTask,
)
from mecfs_bio.build_system.task.magma.magma_gene_set_analysis_task import (
    GENE_SET_ANALYSIS_OUTPUT_STEM_NAME,
)

MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR_EXTRACTED = (
    CopyFileFromDirectoryTask.create_result_table(
        "pgc2022_sch_human_brain_atlas_gene_covar_results_extracted",
        source_directory_task=MAGMA_PGC2022_SCH_HUMAN_BRAIN_ATLAS_GENE_COVAR,
        path_inside_directory=PurePath(
            str(GENE_SET_ANALYSIS_OUTPUT_STEM_NAME + ".gsa.out")
        ),
        extension=".txt",
        read_spec=DataFrameReadSpec(DataFrameWhiteSpaceSepTextFormat(comment_code="#")),
    )
)
