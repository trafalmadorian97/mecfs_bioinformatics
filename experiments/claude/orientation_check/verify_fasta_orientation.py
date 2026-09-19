"""Verify the orientation of the baseline-LF annotations and the Broad UKBB LD
panel against the hg19 FASTA reference.

Context: the SUSIE-polyfun cleanup wants to drop the runtime allele-orientation
reconcilers (the unordered allele key and HarmonizeGWASWithReferenceViaAlleles)
and instead join everything on (chrom, pos, ea, nea), trusting that every input
is already oriented to a single FASTA reference (NEA == reference allele). That
is only safe if the annotations and the LD panel really are reference-oriented.
This script measures that, genome-wide for the annotations and on the LD label
files that happen to be present locally.

Two checks, following the cheap-transitivity plan:

  A. annotations vs FASTA (genome-wide). For every annotation row, is A1 or A2
     the reference base at that position? A consistent column (e.g. A2 == REF for
     ~100% of rows) means the annotation matrix has a fixed reference orientation.

  B. local LD labels vs FASTA and vs annotations. The full LD matrices are not
     stored locally, but their variant label files (rsid, chromosome, position,
     allele1, allele2) are, for whichever intervals have been fine-mapped. For
     those, we (1) check allele1/allele2 against the FASTA directly, and (2) join
     to the annotation matrix on (CHR, POS) to confirm the LD variants are the
     same variants with the same allele set and column order. Together these
     validate the assumption that the LD panel shares the annotations' (and hence
     the FASTA's) orientation, rather than assuming the transitivity.

The generator renames the panel columns allele1 -> non-effect allele, allele2 ->
effect allele (polyfun issue #208), so "panel is reference-oriented with
NEA == REF" means allele1 == REF genome-wide.

Run: pixi r python experiments/claude/orientation_check/verify_fasta_orientation.py
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import polars as pl
import yaml

from mecfs_bio.assets.reference_data.genome_sequence.ucsc_hg19_fasta import (
    UCSC_HG19_INDEXED_FASTA,
)
from mecfs_bio.assets.gwas.ukbb_ppp.ppp_database.hapmap3.hapmap_3_ppp_index import (
    HAPMAP_3_PPP_DATABASE_INDEX,
)
from mecfs_bio.assets.reference_data.polyfun.annotations.baseline_lf_annotations import (
    BASELINE_LF_ANNOTATION_MATRIX,
)
from mecfs_bio.assets.reference_data.polyfun.precomputed_prior.polyfun_precomputed_prior import (
    COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS,
)
from mecfs_bio.build_system.rebuilder.metadata_to_path.remapping_meta_to_path import (
    PathRemapRule,
    RemappingMetaToPath,
)
from mecfs_bio.build_system.task.genome_reference_harmonization.fasta import (
    IndexedFasta,
    reference_matches,
)

_CONFIG_PATH = Path("default_runner_config.yaml")
_LOG_PATH = Path("experiments/claude/orientation_check/verify_fasta_orientation.log")

# Annotation columns (baselineLF2.2.UKB.<chr>.annot.parquet -> matrix).
_ANNOT_CHR = "CHR"
_ANNOT_BP = "BP"
_ANNOT_A1 = "A1"
_ANNOT_A2 = "A2"

# LD label columns (Broad UKBB <interval>.gz).
_LD_CHR = "chromosome"
_LD_POS = "position"
_LD_A1 = "allele1"  # renamed to non-effect allele downstream
_LD_A2 = "allele2"  # renamed to effect allele downstream


class _Tee:
    """Write to stdout and a logfile at once, so the run is reproducible from the log."""

    def __init__(self, log_path: Path) -> None:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        self._log = open(log_path, "w")

    def __call__(self, *parts: object) -> None:
        line = " ".join(str(p) for p in parts)
        print(line)
        self._log.write(line + "\n")
        self._log.flush()

    def close(self) -> None:
        self._log.close()


def _build_resolver() -> RemappingMetaToPath:
    cfg = yaml.safe_load(_CONFIG_PATH.read_text())
    return RemappingMetaToPath(
        default_root=Path(cfg["asset_root"]),
        rules=PathRemapRule.tuple_from_config(cfg.get("path_remap", [])),
    )


def _is_snv(allele_col: str) -> pl.Expr:
    return pl.col(allele_col).str.len_bytes() == 1


def _match_fraction(
    fasta: IndexedFasta, chrom: int, positions: np.ndarray, alleles: pl.Series
) -> np.ndarray:
    """Boolean array: does each allele equal the reference at its position."""
    return reference_matches(
        fasta=fasta, chrom=chrom, positions=positions, alleles=alleles
    )


def check_a1a2_orientation(
    emit: _Tee,
    fasta: IndexedFasta,
    path: Path,
    *,
    label: str,
    chr_col: str,
    pos_col: str,
    a1_col: str,
    a2_col: str,
) -> None:
    """Genome-wide: for a hg19 A1/A2 frame, what fraction has A1 (resp. A2) equal to
    the reference base at its position, and how many rows are indels."""
    emit(f"\n=== {label} vs FASTA (genome-wide) ===")
    chroms = (
        pl.scan_parquet(path)
        .select(chr_col)
        .unique()
        .collect()
        .to_series()
        .sort()
        .to_list()
    )
    emit(f"chromosomes: {chroms}")

    tot = a1_ref = a2_ref = neither = indel_rows = 0
    for chrom in chroms:
        frame = (
            pl.scan_parquet(path)
            .filter(pl.col(chr_col) == chrom)
            .select(
                pl.col(chr_col),
                pl.col(pos_col).cast(pl.Int64),
                pl.col(a1_col).cast(pl.String),
                pl.col(a2_col).cast(pl.String),
            )
            .collect()
        )
        pos = frame[pos_col].to_numpy()
        a1_hit = _match_fraction(fasta, int(chrom), pos, frame[a1_col])
        a2_hit = _match_fraction(fasta, int(chrom), pos, frame[a2_col])
        n = frame.height
        n_a1 = int(a1_hit.sum())
        n_a2 = int(a2_hit.sum())
        n_neither = int((~a1_hit & ~a2_hit).sum())
        n_indel = int(
            frame.select((~_is_snv(a1_col) | ~_is_snv(a2_col)).sum()).item()
        )
        tot += n
        a1_ref += n_a1
        a2_ref += n_a2
        neither += n_neither
        indel_rows += n_indel
        emit(
            f"  chr{chrom:>2}: n={n:>9,}  {a1_col}==REF={n_a1 / n:6.2%}  "
            f"{a2_col}==REF={n_a2 / n:6.2%}  neither={n_neither / n:6.2%}  "
            f"indels={n_indel:,}"
        )

    emit(
        f"TOTAL {label}: n={tot:,}\n"
        f"  {a1_col}==REF : {a1_ref / tot:.4%}\n"
        f"  {a2_col}==REF : {a2_ref / tot:.4%}\n"
        f"  neither : {neither / tot:.4%}\n"
        f"  indel rows ({a1_col} or {a2_col} length>1): "
        f"{indel_rows:,} ({indel_rows / tot:.4%})"
    )
    ref_col = a2_col if a2_ref >= a1_ref else a1_col
    emit(
        f"=> reference-oriented column looks like {ref_col} (higher REF match); "
        f"for an exact join treat {ref_col} as NEA. "
        f"(A prefix-match inflates the longer allele's REF rate on indels, so read "
        f"the clean SNV signal: one column ~100%, neither ~0%.)"
    )


def count_ppp_index_indels(emit: _Tee, index_path: Path) -> None:
    """The PPP database infrastructure keys on the unordered {EA, NEA} set, which is
    unsafe for indels (mirrored T/TCA vs TCA/T collide). Count how many index rows
    are indels -- these are exactly what a SNV assertion on the index would reject."""
    emit("\n=== PPP HAPMAP3 INDEX: indel count ===")
    frame = (
        pl.scan_parquet(index_path)
        .select(
            pl.col("EA").cast(pl.String),
            pl.col("NEA").cast(pl.String),
            pl.col("is_strand_ambiguous"),
        )
        .collect()
    )
    n = frame.height
    is_indel = (frame["EA"].str.len_bytes() > 1) | (frame["NEA"].str.len_bytes() > 1)
    n_indel = int(is_indel.sum())
    n_snv = n - n_indel
    emit(
        f"index rows: {n:,}\n"
        f"  SNV rows   : {n_snv:,} ({n_snv / n:.4%})\n"
        f"  indel rows : {n_indel:,} ({n_indel / n:.4%})"
    )
    if n_indel:
        emit("  example indel rows:")
        for row in frame.filter(is_indel).head(5).iter_rows(named=True):
            emit(f"    EA={row['EA']!r} NEA={row['NEA']!r}")


def _load_local_ld_labels(emit: _Tee, resolver: RemappingMetaToPath) -> pl.DataFrame:
    ld_root = resolver.default_root / "reference_data" / "ukbb_reference_ld"
    gz_files = sorted(ld_root.glob("*/raw/*.gz"))
    emit(f"local LD label files: {len(gz_files)}")
    frames = []
    for gz in gz_files:
        with gzip.open(gz, "rt") as fh:
            frame = pl.read_csv(fh.read().encode(), separator="\t")
        emit(f"  {gz.parent.parent.name}: {frame.height:,} variants")
        frames.append(
            frame.select(
                pl.col(_LD_CHR).cast(pl.Int64),
                pl.col(_LD_POS).cast(pl.Int64),
                pl.col(_LD_A1).cast(pl.String),
                pl.col(_LD_A2).cast(pl.String),
            )
        )
    combined = pl.concat(frames, how="vertical").unique()
    emit(f"combined unique LD variants: {combined.height:,}")
    return combined


def check_ld_vs_fasta_and_annotations(
    emit: _Tee, fasta: IndexedFasta, annot_path: Path, resolver: RemappingMetaToPath
) -> None:
    emit("\n=== B. LOCAL LD LABELS vs FASTA and vs ANNOTATIONS ===")
    ld = _load_local_ld_labels(emit, resolver)

    # B1: LD alleles vs FASTA directly (ground truth for the panel's orientation).
    emit("\n-- B1. LD labels vs FASTA --")
    tot = a1_ref = a2_ref = neither = 0
    for chrom in sorted(ld[_LD_CHR].unique().to_list()):
        sub = ld.filter(pl.col(_LD_CHR) == chrom)
        pos = sub[_LD_POS].to_numpy()
        a1_hit = _match_fraction(fasta, int(chrom), pos, sub[_LD_A1])
        a2_hit = _match_fraction(fasta, int(chrom), pos, sub[_LD_A2])
        n = sub.height
        tot += n
        a1_ref += int(a1_hit.sum())
        a2_ref += int(a2_hit.sum())
        neither += int((~a1_hit & ~a2_hit).sum())
        emit(
            f"  chr{chrom:>2}: n={n:>8,}  allele1==REF={a1_hit.mean():6.2%}  "
            f"allele2==REF={a2_hit.mean():6.2%}"
        )
    emit(
        f"TOTAL local LD: n={tot:,}\n"
        f"  allele1==REF : {a1_ref / tot:.4%}\n"
        f"  allele2==REF : {a2_ref / tot:.4%}\n"
        f"  neither      : {neither / tot:.4%}\n"
        f"(allele1 is the panel non-effect allele; ~100% allele1==REF means the "
        f"panel is reference-oriented with NEA==REF)"
    )

    # B2: LD variants vs annotation variants (identity + allele-column order).
    emit("\n-- B2. LD labels vs annotations (join on CHR, POS) --")
    ld_chroms = sorted(ld[_LD_CHR].unique().to_list())
    annot = (
        pl.scan_parquet(annot_path)
        .filter(pl.col(_ANNOT_CHR).is_in(ld_chroms))
        .select(_ANNOT_CHR, _ANNOT_BP, _ANNOT_A1, _ANNOT_A2)
        .collect()
    )
    joined = ld.join(
        annot,
        left_on=[_LD_CHR, _LD_POS],
        right_on=[_ANNOT_CHR, _ANNOT_BP],
        how="left",
    )
    covered = joined.filter(pl.col(_ANNOT_A1).is_not_null())
    n_ld = joined.height
    n_cov = covered.height
    same_set = covered.filter(
        pl.min_horizontal(_LD_A1, _LD_A2) == pl.min_horizontal(_ANNOT_A1, _ANNOT_A2)
    ).filter(
        pl.max_horizontal(_LD_A1, _LD_A2) == pl.max_horizontal(_ANNOT_A1, _ANNOT_A2)
    )
    n_same_set = same_set.height
    # Among same-allele-set variants, does the panel column order match A1/A2 order?
    same_order = same_set.filter(
        (pl.col(_LD_A1) == pl.col(_ANNOT_A1)) & (pl.col(_LD_A2) == pl.col(_ANNOT_A2))
    ).height
    emit(
        f"LD variants: {n_ld:,}\n"
        f"  covered by annotations (same CHR,POS): {n_cov:,} ({n_cov / n_ld:.4%})\n"
        f"  of covered, same allele SET {{a1,a2}}=={{A1,A2}}: "
        f"{n_same_set:,} ({n_same_set / max(n_cov, 1):.4%})\n"
        f"  of same-set, same column order (allele1==A1 & allele2==A2): "
        f"{same_order:,} ({same_order / max(n_same_set, 1):.4%})"
    )
    emit(
        "(high coverage + same-set ~100% confirms the LD panel and annotations are "
        "the same variants; the column-order fraction tells you whether the panel's "
        "allele1/allele2 follows the annotations' A1/A2 order or the reverse.)"
    )


def main() -> int:
    emit = _Tee(_LOG_PATH)
    try:
        resolver = _build_resolver()
        fasta_dir = resolver(UCSC_HG19_INDEXED_FASTA.meta)
        annot_path = resolver(BASELINE_LF_ANNOTATION_MATRIX.meta)
        emit(f"FASTA dir   : {fasta_dir}")
        emit(f"annotations : {annot_path}")
        assert fasta_dir.exists(), fasta_dir
        assert annot_path.exists(), annot_path
        fasta = IndexedFasta.open(fasta_dir)
        emit(f"FASTA chromosomes indexed: {sorted(fasta.entries)}")

        check_a1a2_orientation(
            emit,
            fasta,
            annot_path,
            label="A. ANNOTATIONS",
            chr_col=_ANNOT_CHR,
            pos_col=_ANNOT_BP,
            a1_col=_ANNOT_A1,
            a2_col=_ANNOT_A2,
        )
        check_ld_vs_fasta_and_annotations(emit, fasta, annot_path, resolver)

        # C. polyfun prior / snpvar meta (same asset): drives the SUSIE prior join
        # and the ridge-weights snpvar join. hg19, columns CHR/BP/A1/A2.
        prior_path = resolver(COMBINED_POLYFUN_PRECOMPUTED_HERITABILITY_WEIGHTS.meta)
        emit(f"\nprior/snpvar: {prior_path}")
        assert prior_path.exists(), prior_path
        check_a1a2_orientation(
            emit,
            fasta,
            prior_path,
            label="C. POLYFUN PRIOR / SNPVAR META",
            chr_col="CHR",
            pos_col="BP",
            a1_col="A1",
            a2_col="A2",
        )

        # D. PPP HapMap3 index indel count (hg38 primary; orientation not checked here,
        # only the indel count relevant to the unordered-key SNV assertion).
        index_path = resolver(HAPMAP_3_PPP_DATABASE_INDEX.meta)
        emit(f"\nPPP index: {index_path}")
        assert index_path.exists(), index_path
        count_ppp_index_indels(emit, index_path)
        emit("\nDONE")
    finally:
        emit.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
