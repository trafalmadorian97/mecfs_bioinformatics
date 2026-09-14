"""Audit every data context that is keyed with unordered_allele_key for mirrored-allele pairs.

A "collision group" is a set of >1 rows sharing (CHR, POS, unordered allele key). Each group
is classified as:
  - exact_dup: the rows have the same ORDERED alleles (plain duplicate rows);
  - mirrored_snv: swapped alleles, both single-nucleotide (same variant, orientation flip);
  - mirrored_indel: swapped alleles, at least one allele longer than 1 (e.g. T/TCA vs TCA/T),
    which in these sources are genuinely DISTINCT variants.

Run all sections, or a subset by number:
  pixi r python experiments/claude/polyfun_annot_dedup/audit_mirrored_indels.py [1 2 ...] \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/audit_mirrored_indels.log
"""

import sys
from pathlib import Path

import polars as pl

from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

STORE = Path("assets/base_asset_store")
DECODE = STORE / "gwas/ME_CFS/DecodeME"
ANNOT_MEMBERS = (
    STORE
    / "reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annot_parquet_members"
)
ANNOT_MERGED = (
    STORE / "reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annotations.parquet"
)
PRIOR = (
    STORE
    / "reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet"
)
K = ["_chr", "_pos", "_k"]


def normalize(lazy: pl.LazyFrame, chr_col: str, pos_col: str, a: str, b: str) -> pl.LazyFrame:
    return lazy.with_columns(
        pl.col(chr_col).cast(pl.Utf8).str.replace("^chr", "").alias("_chr"),
        pl.col(pos_col).cast(pl.Int64).alias("_pos"),
        pl.col(a).cast(pl.Utf8).alias("_a"),
        pl.col(b).cast(pl.Utf8).alias("_b"),
    ).with_columns(unordered_allele_key("_a", "_b").alias("_k"))


def collision_summary(label: str, lazy: pl.LazyFrame) -> pl.DataFrame:
    """Print the collision classification for one dataset and return its mirrored-indel keys."""
    groups = (
        lazy.group_by(K)
        .agg(
            pl.len().alias("n"),
            pl.struct("_a", "_b").n_unique().alias("n_ordered"),
            (pl.col("_a").str.len_chars().max() > 1).alias("indel_a"),
            (pl.col("_b").str.len_chars().max() > 1).alias("indel_b"),
        )
        .filter(pl.col("n") > 1)
        .collect()
    )
    total = lazy.select(pl.len()).collect().item()
    indel = groups["indel_a"] | groups["indel_b"]
    mirrored = groups["n_ordered"] > 1
    n_indel_rows = (
        lazy.filter((pl.col("_a").str.len_chars() > 1) | (pl.col("_b").str.len_chars() > 1))
        .select(pl.len())
        .collect()
        .item()
    )
    print(
        f"{label}: rows={total:,} indel_rows={n_indel_rows:,} collision_groups={groups.height:,} "
        f"exact_dup={int((~mirrored).sum()):,} mirrored_snv={int((mirrored & ~indel).sum()):,} "
        f"mirrored_indel={int((mirrored & indel).sum()):,}"
    )
    return groups.filter(mirrored & indel).select(K)


def collision_summary_by_chrom(label: str, lazy: pl.LazyFrame, chroms: list[str]) -> pl.DataFrame:
    """collision_summary run one chromosome at a time (bounded memory for genome-wide tables)."""
    keys = [collision_summary(f"{label} chr{c}", lazy.filter(pl.col("_chr") == c)) for c in chroms]
    out = pl.concat(keys)
    print(f"{label}: TOTAL mirrored_indel groups={out.height:,}")
    return out


def show_rows(label: str, lazy: pl.LazyFrame, keys: pl.DataFrame, extra: list[str]) -> None:
    """Print the full rows behind the given keys."""
    rows = lazy.join(keys.lazy(), on=K, how="inner").sort(K).collect()
    if rows.height:
        print(f"-- {label}: rows on these keys --")
        with pl.Config(tbl_cols=-1, tbl_width_chars=250, tbl_rows=50):
            print(rows.select(*K, "_a", "_b", *extra))


def polyfun_bad_keys() -> pl.DataFrame:
    prior = normalize(pl.scan_parquet(PRIOR), "CHR", "BP", "A1", "A2")
    return (
        prior.group_by(K)
        .agg(pl.struct("_a", "_b").n_unique().alias("n"))
        .filter(pl.col("n") > 1)
        .select(K)
        .collect()
    )


def section_1_2() -> None:
    print("=== 1. PolyFun baseline-LF annotation (BuildBaselineLFAnnotationParquetTask input) ===")
    annot_members = pl.concat(
        [
            normalize(pl.scan_parquet(p), "CHR", "BP", "A1", "A2")
            for p in sorted(ANNOT_MEMBERS.glob("*.annot.parquet"))
        ]
    )
    annot_keys = collision_summary("annot members (all chr)", annot_members)
    collision_summary(
        "annot merged (local, built with old keep=first)",
        normalize(pl.scan_parquet(ANNOT_MERGED), "CHR", "BP", "A1", "A2"),
    )
    print("\n=== 2. PolyFun precomputed prior (ridge snpvar meta + SUSIE prior) ===")
    prior_keys = collision_summary(
        "precomputed prior", normalize(pl.scan_parquet(PRIOR), "CHR", "BP", "A1", "A2")
    )
    print(
        "prior mirrored-indel keys == annot mirrored-indel keys:",
        prior_keys.sort(K).equals(annot_keys.sort(K)),
    )


def section_3() -> None:
    print("\n=== 3. UKBB Broad LD labels (SUSIE reference / harmonization reference) ===")
    for p in sorted(STORE.glob("reference_data/ukbb_reference_ld/*/processed/*labels_renamed.parquet")):
        lazy = normalize(pl.scan_parquet(p), "CHR", "POS", "EA", "NEA")
        keys = collision_summary(p.parent.parent.name, lazy)
        show_rows(p.parent.parent.name, lazy, keys, ["rsID"])


def section_4() -> None:
    print("\n=== 4. DecodeME sumstats (GWAS side of SUSIE prior join and contrast join) ===")
    lazy = normalize(
        pl.scan_parquet(DECODE / "processed/decode_me_gwas_1_liftover_to_37_parquet_file.parquet"),
        "CHR",
        "POS",
        "EA",
        "NEA",
    )
    keys = collision_summary("DecodeME build-37 liftover (genome-wide)", lazy)
    show_rows("DecodeME build-37 liftover", lazy, keys, ["SNPID", "BETA", "EAF"])
    for p in sorted((DECODE / "processed").glob("*harmonized_with_ref.parquet")):
        lazy = normalize(pl.scan_parquet(p), "CHR", "POS", "EA", "NEA")
        keys = collision_summary(p.stem, lazy)
        show_rows(p.stem, lazy, keys, ["SNPID", "BETA", "EAF"])


def section_5() -> None:
    print("\n=== 5. SUSIE run outputs: duplicates + overlap with PolyFun mirrored-indel keys ===")
    bad_keys = polyfun_bad_keys()
    for run in sorted((DECODE / "analysis").glob("*susie*")):
        fg = run / "filtered_gwas.parquet"
        if not fg.exists():
            continue
        lazy = normalize(pl.scan_parquet(fg), "CHR", "POS", "EA", "NEA")
        dup = lazy.select(pl.len() - pl.struct(K).n_unique()).collect().item()
        overlap = lazy.join(bad_keys.lazy(), on=K, how="inner").select(pl.len()).collect().item()
        c, lo, hi = lazy.select(
            pl.col("_chr").first(), pl.col("_pos").min().alias("lo"), pl.col("_pos").max().alias("hi")
        ).collect().row(0)
        in_window = bad_keys.filter((pl.col("_chr") == c) & pl.col("_pos").is_between(lo, hi)).height
        print(
            f"{run.name}: dup_keys={dup} gwas_rows_on_polyfun_mirrored_keys={overlap} "
            f"polyfun_mirrored_pairs_in_window={in_window}"
        )


def section_6() -> None:
    print("\n=== 6. UKB-PPP (ConstructPppVariantIndexTask / align_protein_to_index) ===")
    template = normalize(
        pl.scan_parquet(
            STORE / "reference_data/ukbb_ppp/sun_et_al_2023/extracted/ukbb_ppp_rabgap1l_sumstats_stacked.parquet"
        ),
        "CHROM",
        "GENPOS",
        "ALLELE1",
        "ALLELE0",
    )
    keys = collision_summary_by_chrom(
        "RABGAP1L template protein (stacked regenie, hg38)", template, [str(c) for c in range(1, 23)]
    )
    show_rows("RABGAP1L template (first 10 keys)", template, keys.head(10), ["ID", "BETA", "A1FREQ"])
    for label, path in [
        ("HapMap3 membership", "reference_data/ppp_variant_index_membership/hapmap3/processed/hapmap3_membership.parquet"),
        ("PPP HapMap3 variant index", "reference_data/ukbb_ppp_variant_index/hapmap_3_membership_list/processed/ppp_variant_index.parquet"),
        ("CSF HapMap3 variant index", "reference_data/csf_pqtl_variant_index/hapmap_3_membership_list/processed/csf_variant_index.parquet"),
    ]:
        collision_summary(label, normalize(pl.scan_parquet(STORE / path), "CHR", "POS", "EA", "NEA"))


def section_7() -> None:
    print("\n=== 7. gwaslab 1kg dbSNP151 table (attach_rsid dbSNP side; hg19 as hg38 proxy) ===")
    dbsnp = pl.scan_csv(
        Path.home() / ".gwaslab/1kg_dbsnp151_hg19_auto.txt.gz",
        separator="\t",
        schema_overrides={"CHR": pl.Utf8},
    )
    # Gzipped CSV cannot be predicate-pushed per chromosome cheaply; sink to parquet once.
    cache = Path("experiments/claude/polyfun_annot_dedup/_dbsnp151_hg19.parquet")
    if not cache.exists():
        dbsnp.sink_parquet(cache)
    collision_summary_by_chrom(
        "1kg dbSNP151 hg19",
        normalize(pl.scan_parquet(cache), "CHR", "POS", "EA", "NEA"),
        [str(c) for c in range(1, 23)],
    )


SECTIONS = {"1": section_1_2, "3": section_3, "4": section_4, "5": section_5, "6": section_6, "7": section_7}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(SECTIONS):
        SECTIONS[name]()
