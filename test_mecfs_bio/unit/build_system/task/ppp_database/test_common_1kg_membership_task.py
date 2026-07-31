import polars as pl

from mecfs_bio.build_system.task.ppp_database.common_1kg_membership_task import (
    attach_rsid,
)
from mecfs_bio.build_system.task.ppp_database.filter_common_1kg_variants_task import (
    ONEKG_AF_COL,
    ONEKG_ALT_COL,
    ONEKG_REF_COL,
)
from mecfs_bio.build_system.task.ppp_database.hapmap3_membership_task import (
    MEMBERSHIP_COLUMNS,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)


def test_attach_rsid():
    filtered = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1],
            GWASLAB_POS_COL: [100, 200],
            ONEKG_REF_COL: ["A", "C"],
            ONEKG_ALT_COL: ["G", "T"],
            ONEKG_AF_COL: [0.2, 0.3],
        }
    )
    # dbSNP row for chr1:100 is in the swapped orientation (EA=G, NEA=A) and must
    # still match; chr1:200 is absent -> null rsID.
    dbsnp = pl.DataFrame(
        {
            GWASLAB_CHROM_COL: [1, 1],
            GWASLAB_POS_COL: [100, 999],
            GWASLAB_EFFECT_ALLELE_COL: ["G", "A"],
            GWASLAB_NON_EFFECT_ALLELE_COL: ["A", "C"],
            GWASLAB_RSID_COL: ["rs1", "rs9"],
        }
    )

    out = attach_rsid(filtered, dbsnp)

    assert out.columns == MEMBERSHIP_COLUMNS
    assert out[GWASLAB_CHROM_COL].to_list() == [1, 1]
    assert out[GWASLAB_POS_COL].to_list() == [100, 200]
    # EA = ALT, NEA = REF.
    assert out[GWASLAB_EFFECT_ALLELE_COL].to_list() == ["G", "T"]
    assert out[GWASLAB_NON_EFFECT_ALLELE_COL].to_list() == ["A", "C"]
    # rsID attached where matched, null (nullable) where absent.
    assert out[GWASLAB_RSID_COL].to_list() == ["rs1", None]
