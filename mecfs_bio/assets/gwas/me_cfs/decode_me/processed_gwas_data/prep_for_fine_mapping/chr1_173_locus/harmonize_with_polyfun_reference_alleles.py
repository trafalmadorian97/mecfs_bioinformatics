from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED,
)
from mecfs_bio.assets.reference_data.ukbb_ld_matrices.from_polyfun.chr1_173000001.chr1_173000001_176000001_labels import (
    CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD,
)
from mecfs_bio.build_system.task.harmonize_gwas_with_reference_table_via_chrom_pos_alleles import (
    HarmonizeGWASWithReferenceViaAlleles,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.pipes.uniquepipe import UniquePipe
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)

DECODE_ME_HARMONIZE_WITH_CHR1_173_000_001_LD_VIA_ALLELES = HarmonizeGWASWithReferenceViaAlleles.create(
    asset_id="decode_me_1_harmonize_with_ld_chr1_173000001_via_alleles",
    gwas_data_task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED.join_task,
    reference_task=CHR1_173000001_17600000_UKBB_LD_LABELS_DOWNLOAD,
    palindrome_strategy="drop",
    gwas_pipe=CompositePipe(
        [
            RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
            UniquePipe(
                by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                ],
                keep="none",
                order_by=[
                    GWASLAB_CHROM_COL,
                    GWASLAB_POS_COL,
                    GWASLAB_EFFECT_ALLELE_COL,
                    GWASLAB_NON_EFFECT_ALLELE_COL,
                    GWASLAB_RSID_COL,
                ],
            ),
        ]
    ),
    ref_pipe=CompositePipe(
        [
            RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
            RenameColPipe(old_name="chromosome", new_name=GWASLAB_CHROM_COL),
            RenameColPipe(old_name="position", new_name=GWASLAB_POS_COL),
            RenameColPipe(
                old_name="allele1",
                new_name=GWASLAB_NON_EFFECT_ALLELE_COL,  # See: https://github.com/omerwe/polyfun/issues/208#issuecomment-2563832487
            ),
            RenameColPipe(old_name="allele2", new_name=GWASLAB_EFFECT_ALLELE_COL),
        ]
    ),
)
