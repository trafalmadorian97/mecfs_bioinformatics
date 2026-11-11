"""
For Western et al. AD-coloc aptamers that were not candidates in our AD x CSF run: is that
because the AD GWAS has no genome-wide-significant lead near their pQTL (min AD p within
+-500 kb), because the pQTL lies in the excluded MHC, or because the aptamer has no hit segment
(p < 5e-8 over HapMap3) near it in our CSF hit table?

Usage: pixi r python experiments/claude/csf_coloc/ad_miss_diagnosis.py
"""

import math
from pathlib import Path

import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.alzheimers.bellenguez_et_al.processed.bellenguez_build_38_harmonized import (
    BELLENGUEZ_BUILD_38_HARMONIZED,
)
from mecfs_bio.assets.gwas.csf_pqtl.csf_database.hapmap3.hapmap3_csf_hit_table import (
    HAPMAP_3_CSF_HIT_TABLE_5E8,
)
from mecfs_bio.constants.genomic_coordinate_constants import extended_mhc_interval
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_CHROM_COL,
    GWASLAB_P_COL,
    GWASLAB_POS_COL,
)
from mecfs_bio.constants.protein_coloc_constants import (
    HIT_CHR_COL,
    HIT_LEAD_MLOG10P_COL,
    HIT_PROTEIN_ID_COL,
    HIT_SEG_END_COL,
    HIT_SEG_START_COL,
)

HERE = Path(__file__).parent
MHC = extended_mhc_interval("38")
MISSES = [
    "10576-7",
    "3166-92",
    "14054-17",
    "7244-16",
    "3642-4",
    "8644-46",
    "8465-52",
    "4556-10",
    "9197-4",
    "5280-68",
    "5430-66",
    "12562-1",
    "11263-57",
    "8840-61",
    "20203-45",
    "8874-53",
    "9451-20",
    "3007-7",
    "10737-96",
    "10851-77",
    "20533-39",
    "8687-26",
    "3066-12",
    "7009-8",
]
HALF = 500_000


def main() -> None:
    reference = pl.read_csv(HERE / "western_ad_proteins.csv").filter(
        pl.col("seq_id").is_in(MISSES), pl.col("pqtl_id").str.contains(":")
    )
    ad = pl.scan_parquet(
        DEFAULT_RUNNER.meta_to_path(BELLENGUEZ_BUILD_38_HARMONIZED.meta)
    )
    hits = pl.read_parquet(DEFAULT_RUNNER.meta_to_path(HAPMAP_3_CSF_HIT_TABLE_5E8.meta))
    rows = []
    for r in reference.unique(["seq_id", "pqtl_id"]).iter_rows(named=True):
        chrom_s, pos_s = r["pqtl_id"].split(":")[:2]
        chrom, pos = int(chrom_s.removeprefix("chr")), int(pos_s)
        ad_min_p = (
            ad.filter(
                pl.col(GWASLAB_CHROM_COL) == chrom,
                pl.col(GWASLAB_POS_COL).is_between(pos - HALF, pos + HALF),
            )
            .select(pl.col(GWASLAB_P_COL).min())
            .collect()
            .item()
        )
        near = hits.filter(
            pl.col(HIT_PROTEIN_ID_COL) == r["seq_id"],
            pl.col(HIT_CHR_COL) == chrom,
            pl.col(HIT_SEG_END_COL) >= pos - HALF,
            pl.col(HIT_SEG_START_COL) <= pos + HALF,
        )
        rows.append(
            {
                "seq_id": r["seq_id"],
                "gene": r["gene"],
                "pqtl": r["pqtl_id"],
                "type": r["pqtl_type"],
                "their_h4": round(r["their_pp_h4"], 3),
                "ad_min_mlog10p_500kb": None
                if ad_min_p is None
                else round(-math.log10(max(ad_min_p, 1e-300)), 1),
                "our_hit_mlog10p_near": near[HIT_LEAD_MLOG10P_COL].max(),
                "in_mhc": chrom == MHC.chrom and MHC.start <= pos <= MHC.end,
            }
        )
    with pl.Config(tbl_rows=60, tbl_width_chars=200):
        print(pl.DataFrame(rows).sort("ad_min_mlog10p_500kb", descending=True))


if __name__ == "__main__":
    main()
