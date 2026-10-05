"""
Size a gnomAD-based allele-frequency panel before committing to a genome-wide pass.

Throwaway probe for the gnomAD AF reference brainstorm. It streams chr21 of a gnomAD genomes
release (v2.1.1 on GRCh37 or v4.1 on GRCh38) over HTTPS through bcftools once, keeping every
record with its FILTER and the per-ancestry AF/AN fields, then measures in polars:

- record counts: total, PASS, SNV/indel, palindromic SNVs;
- per-group max AN, i.e. 2 x the group's sample size, to confirm the group sizes;
- for each candidate threshold tau, the rows and zstd parquet bytes of the panel kept by
  "FILTER == PASS and max over the main groups of AF_g >= tau" (tau = 0 means AF_g > 0
  in at least one main group), with a crude genome-wide extrapolation by VCF byte share.

Extra groups (NFE subgroups, "oth", v4's "ami", "mid", "remaining") are stored as columns
but excluded from the max, because they are subsets of nfe or small and heterogeneous.

Usage:
    pixi r python experiments/claude/gnomad_af_reference/measure_chr21_extract.py RELEASE WORK_DIR \
        2>&1 | tee experiments/claude/gnomad_af_reference/measure_chr21_extract_RELEASE.log

RELEASE is v2 or v4.
"""

import sys
import time
from pathlib import Path

import polars as pl
from attrs import frozen

from mecfs_bio.util.subproc.run_command import execute_command

_BUCKET = "https://gnomad-public-us-east-1.s3.amazonaws.com/release"
THRESHOLDS = [0.0, 1e-4, 1e-3, 1e-2]
_KEY_COLUMNS = ["CHR", "POS", "REF", "ALT"]
_PALINDROMIC_PAIRS = {("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")}


@frozen(slots=True)
class Release:
    vcf_url: str
    chr21_vcf_gib: float
    all_chromosomes_vcf_gib: float
    main_groups: tuple[str, ...]
    extra_groups: tuple[str, ...]

    @property
    def groups(self) -> tuple[str, ...]:
        return self.main_groups + self.extra_groups


RELEASES = {
    "v2": Release(
        vcf_url=f"{_BUCKET}/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.21.vcf.bgz",
        chr21_vcf_gib=6.12,
        # whole-genome file minus exome-intervals file
        all_chromosomes_vcf_gib=460.93 - 9.70,
        main_groups=("afr", "amr", "asj", "eas", "fin", "nfe"),
        extra_groups=("oth", "nfe_nwe", "nfe_seu", "nfe_onf", "nfe_est"),
    ),
    "v4": Release(
        vcf_url=f"{_BUCKET}/4.1/vcf/genomes/gnomad.genomes.v4.1.sites.chr21.vcf.bgz",
        chr21_vcf_gib=7.23,
        all_chromosomes_vcf_gib=524.38,
        main_groups=("afr", "amr", "asj", "eas", "fin", "nfe", "sas"),
        extra_groups=("ami", "mid", "remaining"),
    ),
}


def af_col(group: str) -> str:
    return f"AF_{group}"


def an_col(group: str) -> str:
    return f"AN_{group}"


def stream_records_to_tsv(release: Release, tsv_path: Path) -> float:
    """Stream chr21 into a TSV of keys, FILTER and per-group AF/AN; return elapsed seconds."""
    fields = "".join(
        f"\\t%INFO/{af_col(g)}\\t%INFO/{an_col(g)}" for g in release.groups
    )
    start = time.monotonic()
    execute_command(
        [
            "bcftools",
            "query",
            "-f",
            f"'%CHROM\\t%POS\\t%REF\\t%ALT\\t%FILTER{fields}\\n'",
            release.vcf_url,
            "-o",
            str(tsv_path),
        ]
    )
    return time.monotonic() - start


def scan_records(release: Release, tsv_path: Path) -> pl.LazyFrame:
    schema: dict[str, pl.DataType] = {
        "CHR": pl.String(),
        "POS": pl.Int32(),
        "REF": pl.String(),
        "ALT": pl.String(),
        "FILTER": pl.String(),
    }
    for group in release.groups:
        schema[af_col(group)] = pl.Float32()
        schema[an_col(group)] = pl.Int32()
    return pl.scan_csv(
        tsv_path, separator="\t", has_header=False, null_values=["."], schema=schema
    )


def is_snv() -> pl.Expr:
    return (pl.col("REF").str.len_bytes() == 1) & (pl.col("ALT").str.len_bytes() == 1)


def is_palindromic_snv() -> pl.Expr:
    pair = pl.concat_str("REF", "ALT")
    return is_snv() & pair.is_in(["".join(p) for p in _PALINDROMIC_PAIRS])


def max_main_af(release: Release) -> pl.Expr:
    return pl.max_horizontal(
        [pl.col(af_col(g)).fill_null(0.0) for g in release.main_groups]
    )


def keep_expr(release: Release, tau: float) -> pl.Expr:
    passing = pl.col("FILTER") == "PASS"
    if tau == 0.0:
        return passing & (max_main_af(release) > 0.0)
    return passing & (max_main_af(release) >= tau)


def report_record_counts(records: pl.DataFrame) -> None:
    print("\n== record counts (chr21, all FILTERs)")
    print(
        records.select(
            pl.len().alias("total"),
            (pl.col("FILTER") == "PASS").sum().alias("pass"),
            is_snv().sum().alias("snv"),
            (~is_snv()).sum().alias("non_snv"),
            is_palindromic_snv().sum().alias("palindromic_snv"),
        )
    )
    print("\n== FILTER values")
    print(records.group_by("FILTER").len().sort("len", descending=True))


def report_group_sizes(release: Release, records: pl.DataFrame) -> None:
    print("\n== per-group max AN (approx 2 x samples), PASS records")
    passing = records.filter(pl.col("FILTER") == "PASS")
    with pl.Config(tbl_rows=-1):
        print(
            pl.DataFrame(
                {
                    "group": release.groups,
                    "max_AN": [passing[an_col(g)].max() for g in release.groups],
                    "median_AN": [passing[an_col(g)].median() for g in release.groups],
                }
            ).with_columns((pl.col("max_AN") / 2).alias("approx_samples"))
        )


def panel_columns(release: Release) -> list[str]:
    return _KEY_COLUMNS + [c for g in release.groups for c in (af_col(g), an_col(g))]


def report_thresholds(release: Release, records: pl.DataFrame, out_dir: Path) -> None:
    print("\n== panel size by threshold (PASS and max main-group AF >= tau)")
    scale = release.all_chromosomes_vcf_gib / release.chr21_vcf_gib
    rows = []
    for tau in THRESHOLDS:
        panel = (
            records.filter(keep_expr(release, tau))
            .select(panel_columns(release))
            .sort("POS")
        )
        path = out_dir / f"chr21_tau_{tau:g}.parquet"
        panel.write_parquet(path, compression="zstd")
        size_mib = path.stat().st_size / 2**20
        rows.append(
            {
                "tau": tau,
                "rows": panel.height,
                "snv": panel.select(is_snv().sum()).item(),
                "palindromic_snv": panel.select(is_palindromic_snv().sum()).item(),
                "non_snv": panel.select((~is_snv()).sum()).item(),
                "parquet_mib": round(size_mib, 1),
                "bytes_per_row": round(size_mib * 2**20 / max(panel.height, 1), 1),
                "genome_wide_gib_est": round(size_mib / 1024 * scale, 2),
                "genome_wide_rows_est_m": round(panel.height / 1e6 * scale, 1),
            }
        )
    with pl.Config(tbl_cols=-1, tbl_width_chars=200):
        print(pl.DataFrame(rows))


def main(release: Release, work_dir: Path) -> None:
    work_dir.mkdir(parents=True, exist_ok=True)
    tsv_path = work_dir / "chr21_records.tsv"
    if not tsv_path.exists():
        seconds = stream_records_to_tsv(release, tsv_path)
        genome_hours = (
            seconds * release.all_chromosomes_vcf_gib / release.chr21_vcf_gib / 3600
        )
        print(
            f"streamed chr21 ({release.chr21_vcf_gib} GiB) in {seconds / 60:.1f} min "
            f"({release.chr21_vcf_gib * 1024 / seconds:.1f} MiB/s); genome-wide at this "
            f"rate: {genome_hours:.1f} h"
        )
    records = scan_records(release, tsv_path).collect()
    report_record_counts(records)
    report_group_sizes(release, records)
    report_thresholds(release, records, work_dir)


if __name__ == "__main__":
    assert len(sys.argv) == 3 and sys.argv[1] in RELEASES, (
        f"usage: measure_chr21_extract.py {{{','.join(RELEASES)}}} WORK_DIR"
    )
    main(RELEASES[sys.argv[1]], Path(sys.argv[2]))
