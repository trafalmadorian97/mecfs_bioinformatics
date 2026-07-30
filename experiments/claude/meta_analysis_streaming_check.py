"""
Does the union meta-analysis pipeline stay out-of-core?

FixedEffectsMetaAnalysisTask replaces an iterative chain of pairwise joins with
canonicalise -> vertical concat -> single group_by -> sort -> sink_parquet.
The requirement is that peak memory be O(output size) rather than O(sum of input
sizes), so that a meta analysis whose constituent studies do not collectively fit in
RAM still runs.

To separate the two, three conditions are run, each in its own subprocess so that peak
RSS is not contaminated by a previous condition:

    A  3 studies x M variants, fully overlapping   ->  3M input rows,  M output rows
    B  3 studies x M variants, disjoint            ->  3M input rows, 3M output rows
    C  3 studies x 3M variants, fully overlapping  ->  9M input rows, 3M output rows

If peak memory is driven by the output, A is much smaller than B, and B is close to C.
If it is driven by the input, A is close to B, and both are much smaller than C.

An earlier revision of this script also ran each condition with the final sort stripped
off, to price it: the sort costs 0.2-0.5 GiB of peak RSS and does not change the shape
of the result. The shipped pipeline always sorts, so the toggle is gone.

Run with:
    pixi r python experiments/claude/meta_analysis_streaming_check.py \
        2>&1 | tee experiments/claude/logs/meta_analysis_streaming_check.log
"""

import argparse
import json
import resource
import subprocess
import sys
from pathlib import Path

import narwhals
import polars as pl

from mecfs_bio.build_system.task.fixed_effect_meta_analysis_task import (
    CaseControlSampleInfo,
    canonicalize_source,
    combine_canonicalized_sources,
)

VARIANTS_PER_STUDY = 4_000_000
NUM_STUDIES = 3

CHR_COL = "CHR"
POS_COL = "POS"
EA_COL = "EA"
NEA_COL = "NEA"
BETA_COL = "BETA"
SE_COL = "SE"



def write_study(path: Path, num_variants: int, pos_offset: int, study_index: int):
    """Generate a synthetic study lazily and sink it, so generation itself never
    materialises the study in memory."""
    pos = pl.int_range(
        pos_offset, pos_offset + num_variants, dtype=pl.Int64
    ).alias(POS_COL)
    # A cheap deterministic scramble, so beta/se vary per variant without an RNG buffer.
    scramble = (pl.col(POS_COL) * 2654435761 + study_index * 97) % 1000
    # Half the variants in studies after the first are written in the opposite allele
    # orientation, so the canonicalisation path is exercised.
    reversed_coding = (pl.col(POS_COL) % 2 == 0) & (study_index > 0)
    lf = (
        pl.LazyFrame({"_seed": [0]})
        .select(pos)
        .with_columns(
            (pl.col(POS_COL) % 22 + 1).cast(pl.Int32).alias(CHR_COL),
            pl.when(reversed_coding).then(pl.lit("C")).otherwise(pl.lit("A")).alias(EA_COL),
            pl.when(reversed_coding).then(pl.lit("A")).otherwise(pl.lit("C")).alias(NEA_COL),
            ((scramble - 500) / 1000).cast(pl.Float64).alias(BETA_COL),
            ((scramble + 1) / 10000).cast(pl.Float64).alias(SE_COL),
        )
        .select(CHR_COL, POS_COL, EA_COL, NEA_COL, BETA_COL, SE_COL)
    )
    lf.sink_parquet(path)


def run_pipeline(study_paths: list[Path], out_path: Path):
    """The shipped pipeline, end to end, all lazy.

    This calls the production canonicalisation and combination directly, so the numbers
    describe the task rather than a stand-in for it.
    """
    frames = [
        canonicalize_source(
            narwhals.scan_parquet(path, backend="polars"),
            source_index=i,
            sample_info=CaseControlSampleInfo(cases=1000 * (i + 1), controls=1000 * (i + 1)),
            keep_rsid=False,
        )
        for i, path in enumerate(study_paths)
    ]
    combine_canonicalized_sources(frames, keep_rsid=False).sink_parquet(out_path)


def peak_rss_gib() -> float:
    # ru_maxrss is in kibibytes on Linux.
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / (1024 * 1024)


CONDITIONS = {
    # name: (variants per study, disjoint?)
    "A_overlap_small": (VARIANTS_PER_STUDY, False),
    "B_disjoint_small": (VARIANTS_PER_STUDY, True),
    "C_overlap_large": (VARIANTS_PER_STUDY * NUM_STUDIES, False),
}


def study_paths_for(work_dir: Path, condition: str) -> list[Path]:
    num_variants, disjoint = CONDITIONS[condition]
    return [
        work_dir / f"{condition}_study_{i}.parquet"
        for i in range(NUM_STUDIES)
    ]


def generate(work_dir: Path, condition: str) -> list[Path]:
    num_variants, disjoint = CONDITIONS[condition]
    paths = study_paths_for(work_dir, condition)
    for i, path in enumerate(paths):
        if path.exists():
            continue
        write_study(
            path,
            num_variants=num_variants,
            pos_offset=i * num_variants if disjoint else 0,
            study_index=i,
        )
    return paths


def run_child(work_dir: Path, condition: str):
    paths = generate(work_dir, condition)
    out_path = work_dir / f"{condition}_out.parquet"
    input_bytes = sum(p.stat().st_size for p in paths)
    run_pipeline(paths, out_path)
    output_rows = (
        pl.scan_parquet(out_path).select(pl.len()).collect(engine="streaming").item()
    )
    print(
        "RESULT "
        + json.dumps(
            {
                "condition": condition,
                "input_rows": sum(CONDITIONS[condition][0] for _ in paths),
                "input_bytes": input_bytes,
                "output_rows": output_rows,
                "output_bytes": out_path.stat().st_size,
                "peak_rss_gib": round(peak_rss_gib(), 3),
            }
        )
    )
    out_path.unlink()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--condition", choices=sorted(CONDITIONS))
    args = parser.parse_args()
    args.work_dir.mkdir(parents=True, exist_ok=True)

    if args.condition is not None:
        run_child(args.work_dir, args.condition)
        return

    results = []
    for condition in sorted(CONDITIONS):
        cmd = [
            sys.executable,
            __file__,
            "--work-dir",
            str(args.work_dir),
            "--condition",
            condition,
        ]
        print(f"--- running {condition} ---", flush=True)
        proc = subprocess.run(cmd, capture_output=True, text=True)
        sys.stdout.write(proc.stdout)
        sys.stderr.write(proc.stderr)
        if proc.returncode != 0:
            print(f"FAILED {condition} rc={proc.returncode}")
            continue
        line = next(l for l in proc.stdout.splitlines() if l.startswith("RESULT "))
        results.append(json.loads(line[len("RESULT ") :]))

    print("\n=== summary ===")
    print(
        f"{'condition':<18}{'in rows':>12}{'out rows':>12}{'peak RSS GiB':>14}"
    )
    for r in results:
        print(
            f"{r['condition']:<18}{r['input_rows']:>12,}"
            f"{r['output_rows']:>12,}{r['peak_rss_gib']:>14.2f}"
        )


if __name__ == "__main__":
    main()
