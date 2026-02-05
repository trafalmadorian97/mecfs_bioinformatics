"""
Apply SUSIE to finemap the main DecodeME GWAS hit on chromosome 1.
As a reference, use the UK Biobank LD matrix from the Broad Institute.
"""

from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.prep_for_fine_mapping.chr1_173_locus.harmonize_with_polyfun_reference_alleles import (
    DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES,
)
from mecfs_bio.assets.reference_data.ukbb_ld_matrices.from_polyfun.chr1_173000001.chr1_173000001_176000001_labels import (
    CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD,
)
from mecfs_bio.assets.reference_data.ukbb_ld_matrices.from_polyfun.chr1_173000001.chr1_173000001_176000001_matrix import (
    CHR1_173000001_17600000_UKBB_LD_MATRIX_DOWNLOAD,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.r_tasks.susie_r_finemap_task import (
    BroadInstituteFormatLDMatrix,
    SusieRFinemapTask,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

DECODE_ME_GWAS_1_SUSIE_FINEMAP_CHR1_173_000_001_LOCUS = SusieRFinemapTask.create(
    asset_id="decode_me_gwas_1_37_susie_finemap_chr1_173000001_locus",
    gwas_data_task=DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES,
    ld_labels_task=CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD,
    ld_matrix_source=BroadInstituteFormatLDMatrix(
        CHR1_173000001_17600000_UKBB_LD_MATRIX_DOWNLOAD
    ),
    effective_sample_size=4 / (1 / 15_579 + 1 / 259_909),  # 4/(1/cases + 1/controls)
    ld_labels_pipe=CompositePipe(
        [
            RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
            RenameColPipe(old_name="chromosome", new_name=GWASLAB_CHROM_COL),
            RenameColPipe(old_name="position", new_name=GWASLAB_POS_COL),
            RenameColPipe(
                old_name="allele1",
                new_name=GWASLAB_NON_EFFECT_ALLELE_COL,  # See: https://github.com/omerwe/polyfun/issues/208#issuecomment-2563832487
            ),
            RenameColPipe(old_name="allele2", new_name=GWASLAB_EFFECT_ALLELE_COL),
        ],
    ),
    subsample=None,
)
