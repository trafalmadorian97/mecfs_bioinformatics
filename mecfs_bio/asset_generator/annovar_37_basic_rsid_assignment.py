"""
Asset generator for assigning rsIDs to genome build-37 GWAS datasets.
"""

import narwhals
from attrs import frozen

from mecfs_bio.assets.reference_data.db_snp.db_sn150_build_37_annovar_proc_parquet_rename_unique import (
    PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_DIRECT_DOWNLOAD,
)
from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.reference_data.gnomad.gnomad_allele_frequency_panels import (
    GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
)
from mecfs_bio.build_system.task.base_task import Task
from mecfs_bio.build_system.task.dataframe_output import (
    ParquetOutFormat,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.genome_reference_harmonization_task import (
    GenomeReferenceHarmonizationTask,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.options import (
    GenomeReferenceHarmonizationOptions,
)
from mecfs_bio.build_system.task.gwaslab.gwaslab_sumstats_to_table_task import (
    GwasLabSumstatsToTableTask,
)
from mecfs_bio.build_system.task.join_dataframes_task import JoinDataFramesTask
from mecfs_bio.build_system.task.pipes.cast_pipe import CastPipe
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.data_processing_pipe import DataProcessingPipe
from mecfs_bio.build_system.task.pipes.drop_indels_pipe import DropIndelsPipe
from mecfs_bio.build_system.task.pipes.identity_pipe import IdentityPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.constants.allele_frequency_panel_constants import GnomadV2Group


@frozen(slots=True)
class RSIDAssignmentTaskGroup:
    """
    Collection of tasks used to assign rsIDs by joining with an existing dataframe of SNPs
    """

    pre_harmonization_table_task: Task
    harmonize_task: Task
    join_task: Task


def annovar_37_basic_rsid_assignment(
    sumstats_task: Task,
    base_name: str,
    panel_ancestry: GnomadV2Group,
    use_gwaslab_rsids_convention: bool = False,
    drop_palindromic_ambiguous: bool = True,
    filter_indels_in_harmonized: bool = False,
) -> RSIDAssignmentTaskGroup:
    """
    Asset generator that creates a chain of tasks to assign rsIDs to existing build 37
    sumstats datasets using the annovar dbSNP reference data.

    The gwaslab Sumstats object is dumped to a table and oriented by genome-reference
    harmonization against the UCSC hg19 FASTA and the gnomAD v2.1.1 genomes panel (the
    rehosted copy, so no 14 h build is needed).

    panel_ancestry is the gnomAD group whose allele frequencies resolve palindromes and
    ambiguous indels. It has no default: choose the group that best matches the GWAS
    cohort, for example nfe_nwe for UK Biobank or other north-western European cohorts, and
    nfe for cohorts spread across Europe. The codes are spelled out in
    mecfs_bio/constants/allele_frequency_panel_constants.py.

    Set drop_palindromic_ambiguous to False to keep palindromic SNVs whose strand cannot be
    resolved. Ambiguous indels are never kept on that basis.

    Set filter_indels_in_harmonized to drop indels before harmonization. This suits datasets
    with very long structural-variant alleles that the downstream SNP-based analyses (LDSC
    genetic correlation, MAGMA) do not use.
    """
    pre_harmonization_table_task = GwasLabSumstatsToTableTask.create_from_source_task(
        source_tsk=sumstats_task,
        asset_id=base_name + "_pre_harmonization_dump_to_parquet",
        sub_dir="processed",
        pipe=DropIndelsPipe() if filter_indels_in_harmonized else IdentityPipe(),
    )
    harmonize_task = GenomeReferenceHarmonizationTask.create(
        asset_id=base_name + "_genome_reference_harmonized",
        sumstats_task=pre_harmonization_table_task,
        fasta_task=UCSC_HG19_INDEXED_FASTA,
        panel_task=GNOMAD_V2_1_1_GENOMES_HG19_ALLELE_FREQUENCIES_REHOSTED,
        panel_ancestry=panel_ancestry,
        options=GenomeReferenceHarmonizationOptions(
            keep_unresolved_palindromes=not drop_palindromic_ambiguous
        ),
    )
    out_pipe: DataProcessingPipe
    if use_gwaslab_rsids_convention:
        out_pipe = RenameColPipe(old_name="rsid", new_name="rsID")
    else:
        out_pipe = IdentityPipe()
    join_with_rsid_task = JoinDataFramesTask.create_from_result_df(
        asset_id=base_name + "_assign_rsids_via_dbsnp150",
        result_df_task=harmonize_task,
        reference_df_task=PARQUET_DBSNP150_37_ANNOVAR_PROC_RENAME_UNIQUE_DIRECT_DOWNLOAD,
        left_on=["CHR", "POS", "EA", "NEA"],
        right_on=["int_chrom", "POS", "ALT", "REF"],
        out_format=ParquetOutFormat(),
        how="inner",
        df_1_pipe=CompositePipe(
            [
                CastPipe(
                    target_column="EA",
                    type=narwhals.dtypes.String(),
                    new_col_name="EA",
                ),
                CastPipe(
                    target_column="NEA",
                    type=narwhals.dtypes.String(),
                    new_col_name="NEA",
                ),
            ]
        ),
        backend="ibis",
        out_pipe=out_pipe,
    )
    return RSIDAssignmentTaskGroup(
        pre_harmonization_table_task=pre_harmonization_table_task,
        harmonize_task=harmonize_task,
        join_task=join_with_rsid_task,
    )
