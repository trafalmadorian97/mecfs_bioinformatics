from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.ppp_database.hapmap3_membership_task import (
    MEMBERSHIP_COLUMNS,
    bundled_hapmap3_hg38_snplist_path,
    normalize_hapmap3_snplist,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_NON_EFFECT_ALLELE_COL,
    GWASLAB_POS_COL,
    GWASLAB_RSID_COL,
)


def test_bundled_hapmap3_snplist_present():
    # gwaslab is a pinned dependency, so its bundled snplist must exist.
    assert bundled_hapmap3_hg38_snplist_path().is_file()


def test_normalize_hapmap3_snplist(tmp_path: Path):
    snplist = tmp_path / "hm3.snplist"
    snplist.write_text(
        "rsid\tA1\tA2\t#CHROM\tPOS\nrs1\tA\tG\t1\t1000\nrs2\tC\tT\t2\t2000\n"
    )

    out = normalize_hapmap3_snplist(snplist)

    assert out.columns == MEMBERSHIP_COLUMNS
    assert out[GWASLAB_CHROM_COL].to_list() == [1, 2]
    assert out[GWASLAB_POS_COL].to_list() == [1000, 2000]
    assert out[GWASLAB_EFFECT_ALLELE_COL].to_list() == ["A", "C"]
    assert out[GWASLAB_NON_EFFECT_ALLELE_COL].to_list() == ["G", "T"]
    assert out[GWASLAB_RSID_COL].to_list() == ["rs1", "rs2"]
    assert out.schema[GWASLAB_CHROM_COL] == pl.Int32
    assert out.schema[GWASLAB_POS_COL] == pl.Int32
