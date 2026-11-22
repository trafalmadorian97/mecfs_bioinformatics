"""
Check whether deCODE's variant_info_iceland.txt.gz can supply an effect-allele frequency (EAF)
for a SomaScan protein sumstats file.

The sumstats carry only ImpMAF. variant_info_iceland has Chrom, Pos, A1, ImpFreqA1, ImpInfo and
no second allele. The design (section 4) wants EAF = ImpFreqA1 where A1 == effectAllele. On one
whole protein file, this script measures:

  1. how variant_info is laid out: exact duplicate rows, rows per position, and whether
     ImpFreqA1 sums to 1 over the distinct A1 values at a position (one row per allele);
  2. for each protein row: how many variant_info rows sit at its (Chrom, Pos), and whether A1
     equals effectAllele (exactly once, more than once), only otherAllele, or neither;
  3. whether the frequencies agree on the A1 == effectAllele rows: min(f, 1 - f) vs ImpMAF,
     with f = ImpFreqA1;
  4. points 2 and 3 again on the shared HapMap3 index rows only (the rows the database keeps).

It runs one chromosome at a time to stay inside 16 GB. Both inputs are first converted to
parquet with streaming, in a cache directory.

Usage (from the repo root):
  pixi r python experiments/claude/decode_proteomics_survey/check_variant_info_alleles.py \
      <protein.txt.gz> <variant_info_iceland.txt.gz> <shared_variant_index.parquet> <cache_dir> \
      | tee experiments/claude/decode_proteomics_survey/check_variant_info_alleles.log
"""

import sys
import time
from pathlib import Path

import polars as pl

T0 = time.perf_counter()


def log(msg: str) -> None:
    print(f"[{time.perf_counter() - T0:6.1f}s] {msg}", flush=True)


def to_parquet(src: Path, dst: Path, columns: list[str], overrides: dict) -> Path:  # noqa: ANN001
    if not dst.exists():
        pl.scan_csv(src, separator="\t", schema_overrides=overrides).select(columns).sink_parquet(dst)
    return dst


def classify(joined: pl.DataFrame) -> pl.DataFrame:
    """One row per protein row: counts of variant_info rows at the position and allele matches."""
    per_row = joined.group_by("row_number").agg(
        n_vi=pl.col("A1").is_not_null().sum(),
        eff_match=(pl.col("A1") == pl.col("effectAllele")).sum(),
        oth_match=(pl.col("A1") == pl.col("otherAllele")).sum(),
        is_bang=(pl.col("otherAllele") == "!").first(),
    )
    return per_row.select(
        rows=pl.len(),
        n_vi_0=(pl.col("n_vi") == 0).sum(),
        n_vi_1=(pl.col("n_vi") == 1).sum(),
        n_vi_2=(pl.col("n_vi") == 2).sum(),
        n_vi_3plus=(pl.col("n_vi") >= 3).sum(),
        eff_exactly_one=(pl.col("eff_match") == 1).sum(),
        eff_more_than_one=(pl.col("eff_match") > 1).sum(),
        only_other=((pl.col("eff_match") == 0) & (pl.col("oth_match") > 0)).sum(),
        neither=((pl.col("eff_match") == 0) & (pl.col("oth_match") == 0) & (pl.col("n_vi") > 0)).sum(),
        bang_rows=pl.col("is_bang").sum(),
        bang_eff_exactly_one=(pl.col("is_bang") & (pl.col("eff_match") == 1)).sum(),
    )


def matched_freqs(prot: pl.DataFrame, vi_unique: pl.DataFrame, chrom: str) -> pl.DataFrame:
    """Protein rows joined to deduplicated variant_info on (Pos, A1 == effectAllele)."""
    m = prot.join(vi_unique, left_on=["Pos", "effectAllele"], right_on=["Pos", "A1"], how="inner")
    return m.select(
        "row_number", pl.lit(chrom).alias("Chrom"), "Pos", "effectAllele", "otherAllele",
        "ImpFreqA1", "ImpMAF",
        diff=(pl.min_horizontal(pl.col("ImpFreqA1"), 1 - pl.col("ImpFreqA1")) - pl.col("ImpMAF"))
        .abs().cast(pl.Float32),
    )


def summarize_freqs(m: pl.DataFrame, n_prot: int) -> pl.DataFrame:
    per_row = m.group_by("row_number").agg(n=pl.len(), spread=pl.col("ImpFreqA1").max() - pl.col("ImpFreqA1").min())
    return pl.DataFrame({
        "protein_rows": [n_prot],
        "rows_with_a_match": [per_row.height],
        "rows_with_several_matches": [per_row.filter(pl.col("n") > 1).height],
        "several_matches_freq_spread_over_1e3": [per_row.filter(pl.col("spread") > 1e-3).height],
    })


def layout(vi: pl.DataFrame, vi_unique: pl.DataFrame) -> pl.DataFrame:
    per_pos = vi_unique.group_by("Pos").agg(n=pl.len(), fsum=pl.col("ImpFreqA1").sum(),
                                            dup_a1=pl.col("A1").n_unique() < pl.len())
    return per_pos.select(
        positions=pl.len(),
        raw_rows=pl.lit(vi.height),
        exact_duplicate_rows=pl.lit(vi.height - vi_unique.height),
        rows=pl.col("n").sum(),
        pos_with_1=(pl.col("n") == 1).sum(),
        pos_with_2=(pl.col("n") == 2).sum(),
        pos_with_3plus=(pl.col("n") >= 3).sum(),
        pos_fsum_within_1e3_of_1=((pl.col("fsum") - 1).abs() < 1e-3).sum(),
        pos_with_repeated_a1_after_dedup=pl.col("dup_a1").sum(),
    )


def report_counts(label: str, df: pl.DataFrame, denom: str) -> None:
    tot = df.sum().row(0, named=True)
    n = tot[denom]
    for k, v in tot.items():
        log(f"[{label}] {k}: {v:,}" + ("" if k == denom else f" ({100 * v / n:.4f}% of {denom})"))


def report_freqs(label: str, m: pl.DataFrame) -> None:
    log(f"[{label}] frequency agreement on A1 == effectAllele rows: {freq_summary(m)}")
    worst = m.top_k(8, by="diff").drop("row_number")
    log(f"[{label}] largest |min(f,1-f) - ImpMAF|:\n{worst}")


def freq_summary(m: pl.DataFrame) -> dict:  # noqa: ANN001
    return m.select(
        n=pl.len(),
        max_abs_diff=pl.col("diff").max(),
        p99=pl.col("diff").quantile(0.99),
        p999=pl.col("diff").quantile(0.999),
        over_1e3=(pl.col("diff") > 1e-3).sum(),
        over_1e2=(pl.col("diff") > 1e-2).sum(),
        effect_allele_is_minor=(pl.col("ImpFreqA1") <= 0.5).mean(),
    ).row(0, named=True)


def report_freqs_light(label: str, m: pl.DataFrame) -> None:
    log(f"[{label}] frequency agreement on A1 == effectAllele rows: {freq_summary(m)}")


def main() -> None:
    prot_path, vi_path, index_path, cache = (Path(a) for a in sys.argv[1:5])
    cache.mkdir(parents=True, exist_ok=True)
    prot_pq = to_parquet(prot_path, cache / "protein.parquet",
                         ["Chrom", "Pos", "Name", "effectAllele", "otherAllele", "ImpMAF"],
                         {"Pos": pl.Int64, "ImpMAF": pl.Float64})
    vi_pq = to_parquet(vi_path, cache / "variant_info.parquet", ["Chrom", "Pos", "A1", "ImpFreqA1"],
                       {"Pos": pl.Int64, "A1": pl.String, "ImpFreqA1": pl.Float64})
    log("parquet conversions ready")

    idx = pl.read_parquet(index_path, columns=["row_number", "Name"]).with_columns(
        pl.col("row_number").cast(pl.UInt32))
    prot_all = pl.scan_parquet(prot_pq).with_row_index("row_number")
    name_ok = (prot_all.join(idx.lazy(), on="row_number").filter(pl.col("Name") == pl.col("Name_right"))
               .select(pl.len()).collect().item())
    log(f"index row_number convention: {name_ok:,} of {idx.height:,} names match")
    assert name_ok == idx.height
    idx_rows = idx["row_number"].implode()

    chroms = pl.scan_parquet(prot_pq).select(pl.col("Chrom").unique(maintain_order=True)).collect()["Chrom"]
    lay, c_all, c_idx, s_all, s_idx, m_all, m_idx = [], [], [], [], [], [], []
    for chrom in chroms:
        prot = prot_all.filter(pl.col("Chrom") == chrom).drop("Name", "Chrom").collect()
        vi = pl.scan_parquet(vi_pq).filter(pl.col("Chrom") == chrom).drop("Chrom").collect()
        vi_unique = vi.unique()
        lay.append(layout(vi, vi_unique))
        joined = prot.join(vi_unique, on="Pos", how="left")
        in_idx = joined.filter(pl.col("row_number").is_in(idx_rows))
        c_all.append(classify(joined))
        c_idx.append(classify(in_idx))
        del joined, in_idx
        m = matched_freqs(prot, vi_unique, chrom)
        m_i = m.filter(pl.col("row_number").is_in(idx_rows))
        n_idx = prot.filter(pl.col("row_number").is_in(idx_rows)).height
        s_all.append(summarize_freqs(m, prot.height))
        s_idx.append(summarize_freqs(m_i, n_idx))
        m_all.append(m.select("diff", "ImpFreqA1"))
        m_idx.append(m_i)
        log(f"{chrom}: protein {prot.height:,} rows, variant_info {vi.height:,} rows")

    report_counts("variant_info layout", pl.concat(lay), "positions")
    report_counts("all protein rows", pl.concat(c_all), "rows")
    report_counts("all protein rows, join on (Pos, A1 == effectAllele)", pl.concat(s_all), "protein_rows")
    report_freqs_light("all protein rows", pl.concat(m_all))
    report_counts("HapMap3 index rows", pl.concat(c_idx), "rows")
    report_counts("HapMap3 index rows, join on (Pos, A1 == effectAllele)", pl.concat(s_idx), "protein_rows")
    report_freqs("HapMap3 index rows", pl.concat(m_idx))


if __name__ == "__main__":
    main()
