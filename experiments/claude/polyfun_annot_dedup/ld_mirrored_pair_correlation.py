"""For each mirrored indel pair in the Broad UKBB LD label files, report the LD correlation between
the two entries (r ~ -1 => one variant listed in both orientations; |r| clearly < 1 => two distinct
variants), plus the PolyFun prior rows and DecodeME build-37 rows at that position."""
from pathlib import Path

import numpy as np
import polars as pl
import scipy.sparse

from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

STORE = Path("assets/base_asset_store")
PRIOR = STORE / "reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet"
DECODE37 = STORE / "gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_liftover_to_37_parquet_file.parquet"

for raw_dir in sorted(STORE.glob("reference_data/ukbb_reference_ld/*/raw")):
    gz = next(raw_dir.glob("*.gz"))
    npz = next(raw_dir.glob("*.npz"))
    labels = pl.read_csv(gz, separator="\t").with_row_index("idx")
    labels = labels.with_columns(unordered_allele_key("allele1", "allele2").alias("k"))
    pairs = labels.filter(pl.len().over(["chromosome", "position", "k"]) > 1)
    if pairs.height == 0:
        continue
    m = scipy.sparse.load_npz(npz).tocsr()
    for (pos, k), grp in pairs.group_by(["position", "k"], maintain_order=True):
        i, j = sorted(grp["idx"].to_list())
        r = m[i, j] + m[j, i]  # upper-triangular storage
        print(f"\n{raw_dir.parent.name} pos={pos} key={k} rsid={grp['rsid'][0]} idx=({i},{j}) r={r:.4f}")
        print(grp.select("idx", "rsid", "allele1", "allele2"))
        chrom = int(grp["chromosome"][0])
        print(pl.scan_parquet(PRIOR).filter((pl.col("CHR") == chrom) & (pl.col("BP") == pos)).collect())
        print(pl.scan_parquet(DECODE37).filter((pl.col("CHR") == chrom) & (pl.col("POS") == pos)).select("SNPID", "EA", "NEA", "EAF", "BETA", "SE").collect())
