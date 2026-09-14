"""Check allele columns against the reference FASTA for every source that carries mirrored indel
pairs (same CHR/POS, swapped alleles, e.g. C/CTAAA and CTAAA/C).

For each source two checks are run.

1. Convention (chr22 sample, SNVs and indels separately): for each row, does column A and/or
   column B match the reference genome starting at POS? This identifies which column (if any)
   is consistently the reference allele.

2. Mirrored pairs (genome-wide): for each pair of rows sharing (CHR, POS, unordered allele key)
   with swapped alleles, check
     - whether the SHORT and the LONG allele each match the genome at POS
       (long matches => a deletion REF=long and an insertion REF=short are both valid records);
     - in how many of the pair's rows the source's candidate reference column matches the genome.
   A pair where both rows are reference-consistent is an insertion + deletion at one position,
   separable by an order-aware (REF, ALT) key. A pair where only one row is reference-consistent
   means one row carries a non-reference ordering, so an order-aware key would be unreliable.

Run:
  pixi r python experiments/claude/polyfun_annot_dedup/fasta_ref_allele_check.py \
    2>&1 | tee experiments/claude/polyfun_annot_dedup/fasta_ref_allele_check.log
"""

from collections import Counter
from collections.abc import Callable
from pathlib import Path

import polars as pl
import pysam
from attrs import frozen

from mecfs_bio.build_system.task.ppp_database.allele_key import unordered_allele_key

STORE = Path("assets/base_asset_store")
GWASLAB = Path.home() / ".gwaslab"
FASTA = {"hg19": GWASLAB / "hg19.fa", "hg38": GWASLAB / "hg38.fa"}
CHROMS = [str(c) for c in range(1, 23)]
SAMPLE_CHROM = "22"
SAMPLE_CAP = 200_000
K = ["chr", "pos", "k"]


@frozen
class Source:
    """A table normalized to columns chr (str), pos (int), a, b, where a is the column we
    believe to be the reference allele. load(chrom) returns that chromosome's rows."""

    name: str
    build: str
    col_a: str
    col_b: str
    load: Callable[[str], pl.DataFrame]


def _norm(lazy: pl.LazyFrame, chr_col: str, pos_col: str, a: str, b: str) -> pl.LazyFrame:
    return lazy.select(
        pl.col(chr_col).cast(pl.Utf8).str.replace("^chr", "").alias("chr"),
        pl.col(pos_col).cast(pl.Int64).alias("pos"),
        pl.col(a).cast(pl.Utf8).str.to_uppercase().alias("a"),
        pl.col(b).cast(pl.Utf8).str.to_uppercase().alias("b"),
    )


def parquet_source(name: str, build: str, path: Path, chr_col: str, pos_col: str, a: str, b: str) -> Source:
    def load(chrom: str) -> pl.DataFrame:
        return (
            _norm(pl.scan_parquet(path), chr_col, pos_col, a, b)
            .filter(pl.col("chr") == chrom)
            .collect()
        )

    return Source(name, build, a, b, load)


def eager_source(name: str, build: str, frame: pl.DataFrame, a: str, b: str) -> Source:
    return Source(name, build, a, b, lambda chrom: frame.filter(pl.col("chr") == chrom))


class Genome:
    def __init__(self, path: Path) -> None:
        self._fasta = pysam.FastaFile(str(path))
        self._chrom: str | None = None
        self._seq = ""

    def seq(self, chrom: str) -> str:
        if chrom != self._chrom:
            self._seq = self._fasta.fetch(f"chr{chrom}").upper()
            self._chrom = chrom
        return self._seq


def _matches(seq: str, pos: int, allele: str) -> bool:
    return seq[pos - 1 : pos - 1 + len(allele)] == allele


def convention(src: Source, genome: Genome) -> None:
    df = src.load(SAMPLE_CHROM)
    if df.height > SAMPLE_CAP:
        df = df.sample(SAMPLE_CAP, seed=0)
    seq = genome.seq(SAMPLE_CHROM)
    counts: dict[str, Counter[str]] = {"SNV": Counter(), "indel": Counter()}
    for pos, a, b in df.select("pos", "a", "b").iter_rows():
        kind = "SNV" if len(a) == 1 and len(b) == 1 else "indel"
        ma, mb = _matches(seq, pos, a), _matches(seq, pos, b)
        label = "both" if ma and mb else f"{src.col_a}_only" if ma else f"{src.col_b}_only" if mb else "neither"
        counts[kind][label] += 1
    for kind, c in counts.items():
        total = sum(c.values())
        if total == 0:
            continue
        parts = ", ".join(f"{k}={v:,} ({100 * v / total:.2f}%)" for k, v in c.most_common())
        print(f"  convention chr{SAMPLE_CHROM} {kind} (n={total:,}): {parts}")


def mirrored_pairs(src: Source, genome: Genome) -> None:
    pair_classes: Counter[str] = Counter()
    long_matches: Counter[str] = Counter()
    examples: list[tuple] = []
    for chrom in CHROMS:
        df = src.load(chrom).with_columns(unordered_allele_key("a", "b").alias("k"))
        pairs = (
            df.filter(
                (pl.col("a").str.len_chars() > 1) | (pl.col("b").str.len_chars() > 1)
            )
            .filter(pl.struct("a", "b").n_unique().over(K) > 1)
            .unique(subset=["chr", "pos", "a", "b"])
        )
        if pairs.height == 0:
            continue
        seq = genome.seq(chrom)
        for (_, pos, _k), grp in pairs.group_by(K, maintain_order=True):
            rows = grp.select("a", "b").rows()
            short, long_ = sorted({rows[0][0], rows[0][1]}, key=len)
            ms, ml = _matches(seq, pos, short), _matches(seq, pos, long_)
            long_matches[f"short={'Y' if ms else 'N'} long={'Y' if ml else 'N'}"] += 1
            n_ref_ok = sum(_matches(seq, pos, a) for a, _ in rows)
            pair_classes[f"{n_ref_ok}/{len(rows)} rows have {src.col_a}==genome"] += 1
            if n_ref_ok < len(rows) and len(examples) < 5:
                examples.append((chrom, pos, rows, seq[pos - 1 : pos - 1 + len(long_) + 5]))
    total = sum(long_matches.values())
    print(f"  mirrored indel pairs genome-wide: {total:,}")
    for k, v in long_matches.most_common():
        print(f"    genome match {k}: {v:,}")
    for k, v in pair_classes.most_common():
        print(f"    {k}: {v:,}")
    for chrom, pos, rows, context in examples:
        print(f"    example not fully ref-consistent: chr{chrom}:{pos} rows(a,b)={rows} genome={context}")


def decode_raw() -> pl.DataFrame:
    """DecodeME raw regenie is gzipped text; read the needed columns once."""
    path = STORE / "gwas/ME_CFS/DecodeME/raw/DecodeME Summary Statistics/gwas_1.regenie.gz"
    return _norm(
        pl.scan_csv(path, separator=" ", schema_overrides={"CHROM": pl.Utf8}),
        "CHROM",
        "GENPOS",
        "ALLELE0",
        "ALLELE1",
    ).collect()


def ld_labels() -> pl.DataFrame:
    """Union of all local Broad UKBB LD label windows (hg19, allele1/allele2)."""
    frames = [
        _norm(pl.scan_csv(p, separator="\t", schema_overrides={"chromosome": pl.Utf8}), "chromosome", "position", "allele1", "allele2")
        for p in sorted(STORE.glob("reference_data/ukbb_reference_ld/*/raw/*.gz"))
    ]
    return pl.concat([f.collect() for f in frames]).unique()


def decode_harmonized() -> pl.DataFrame:
    """Union of the current (palindromes_keep) DecodeME harmonized locus files (hg19, NEA/EA)."""
    paths = sorted((STORE / "gwas/ME_CFS/DecodeME/processed").glob("decode_me_polyfun_explain*harmonized_with_ref.parquet"))
    return pl.concat([_norm(pl.scan_parquet(p), "CHR", "POS", "NEA", "EA").collect() for p in paths]).unique()


def ppp_hg19() -> Callable[[str], pl.DataFrame]:
    path = STORE / "reference_data/ukbb_ppp/sun_et_al_2023/extracted/ukbb_ppp_rabgap1l_sumstats_stacked.parquet"

    def load(chrom: str) -> pl.DataFrame:
        return (
            pl.scan_parquet(path)
            .filter(pl.col("CHROM") == int(chrom))
            .select(
                pl.lit(chrom).alias("chr"),
                pl.col("ID").str.split(":").list.get(1).cast(pl.Int64).alias("pos"),
                pl.col("ALLELE0").str.to_uppercase().alias("a"),
                pl.col("ALLELE1").str.to_uppercase().alias("b"),
            )
            .collect()
        )

    return load


def main() -> None:
    genomes = {build: Genome(path) for build, path in FASTA.items()}
    ppp_path = STORE / "reference_data/ukbb_ppp/sun_et_al_2023/extracted/ukbb_ppp_rabgap1l_sumstats_stacked.parquet"
    annot_dir = STORE / "reference_data/polyfun/annotations/raw/baseline_lf_2.2_ukb_annot_parquet_members"

    def annot_load(chrom: str) -> pl.DataFrame:
        return _norm(pl.scan_parquet(annot_dir / f"baselineLF2.2.UKB.{chrom}.annot.parquet"), "CHR", "BP", "A1", "A2").collect()

    lazy_sources: list[Callable[[], Source]] = [
        lambda: eager_source("DecodeME raw regenie (hg38)", "hg38", decode_raw(), "ALLELE0", "ALLELE1"),
        lambda: parquet_source(
            "DecodeME build-37 liftover (hg19)", "hg19",
            STORE / "gwas/ME_CFS/DecodeME/processed/decode_me_gwas_1_liftover_to_37_parquet_file.parquet",
            "CHR", "POS", "NEA", "EA",
        ),
        lambda: eager_source("DecodeME harmonized polyfun-explain loci (hg19)", "hg19", decode_harmonized(), "NEA", "EA"),
        lambda: parquet_source("UKB-PPP RABGAP1L regenie (hg38 GENPOS)", "hg38", ppp_path, "CHROM", "GENPOS", "ALLELE0", "ALLELE1"),
        lambda: Source("UKB-PPP RABGAP1L regenie (hg19 pos from ID)", "hg19", "ALLELE0", "ALLELE1", ppp_hg19()),
        lambda: parquet_source(
            "PolyFun precomputed prior (hg19)", "hg19",
            STORE / "reference_data/polyfun/precomputed_prior/raw/polyfun_precomputed_heritability_weight_concat.parquet",
            "CHR", "BP", "A1", "A2",
        ),
        lambda: Source("PolyFun baseline-LF annotations (hg19)", "hg19", "A1", "A2", annot_load),
        lambda: eager_source("UKBB Broad LD labels, all local windows (hg19)", "hg19", ld_labels(), "allele1", "allele2"),
    ]
    for make in lazy_sources:
        src = make()
        print(f"\n=== {src.name}: candidate reference column = {src.col_a}, other = {src.col_b} ===", flush=True)
        genome = genomes[src.build]
        convention(src, genome)
        mirrored_pairs(src, genome)


if __name__ == "__main__":
    main()
