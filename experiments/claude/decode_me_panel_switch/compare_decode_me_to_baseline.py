"""
Compare the rebuilt DecodeME outputs (gnomAD v2.1.1 nfe_nwe panel) with the baseline that
record_decode_me_baseline took before the switch (1000 Genomes EUR panel).

Two layers:

- Every baseline asset: are its files byte-identical (md5) to the rebuilt ones? Summarized by
  task type.
- Changed assets of the key analysis types: a content comparison keyed on the natural row
  identifier (variant, cell type, exposure, gene set, ...): rows present on one side only,
  the largest absolute change per numeric column, and significance calls that flipped.
  Fine-mapping additionally compares credible-set membership.

Run after rebuild_decode_me.

Usage:
    pixi r python -m experiments.claude.decode_me_panel_switch.compare_decode_me_to_baseline \
        2>&1 | tee experiments/claude/decode_me_panel_switch/compare_decode_me_to_baseline.log
"""

import collections
import json
from collections.abc import Callable, Sequence
from pathlib import Path

import polars as pl
import yaml
from attrs import frozen

from experiments.claude.decode_me_panel_switch.record_decode_me_baseline import (
    BASELINE_DIR,
    MANIFEST_NAME,
    files_under,
    md5_of,
)
from mecfs_bio.constants.gwaslab_constants import (
    GWASLAB_BETA_COL,
    GWASLAB_EFFECT_ALLELE_COL,
    GWASLAB_EFFECT_ALLELE_FREQ_COL,
    GWASLAB_SNPID_COL,
)

SIGNIFICANCE = 0.05
CORRECTED_P_COL = "_Corrected P Value_"
PIP_ID_COL = "__index_level_0__"
CS_VARIANT_COLS = ["CHR", "POS", "EA", "NEA"]


@frozen(slots=True)
class Record:
    asset_id: str
    task_type: str
    path: Path
    copied: bool
    file_md5: dict[str, str]

    @property
    def new_root(self) -> Path:
        return self.path if self.path.is_dir() else self.path.parent

    @property
    def old_root(self) -> Path:
        return BASELINE_DIR / "outputs" / self.asset_id

    def old(self, name: str) -> Path:
        return self.old_root / name

    def new(self, name: str) -> Path:
        return self.new_root / name


def load_records() -> list[Record]:
    raw = json.loads((BASELINE_DIR / MANIFEST_NAME).read_text())
    return [
        Record(
            asset_id=r["asset_id"],
            task_type=r["task_type"],
            path=Path(r["path"]),
            copied=r["copied"],
            file_md5=r["file_md5"],
        )
        for r in raw
        if r["exists"]
    ]


def new_md5(record: Record) -> dict[str, str] | None:
    if not record.path.exists():
        return None
    root = record.new_root
    return {str(f.relative_to(root)): md5_of(f) for f in files_under(record.path)}


def report_md5_status(records: Sequence[Record]) -> set[str]:
    """Print identical/changed/missing counts per task type; return the changed asset ids."""
    status: dict[str, str] = {}
    changed_files: dict[str, list[str]] = {}
    for record in records:
        current = new_md5(record)
        if current is None:
            status[record.asset_id] = "missing"
        elif current == record.file_md5:
            status[record.asset_id] = "identical"
        else:
            status[record.asset_id] = "changed"
            changed_files[record.asset_id] = sorted(
                name
                for name in set(current) | set(record.file_md5)
                if current.get(name) != record.file_md5.get(name)
            )
    by_type: dict[str, collections.Counter[str]] = collections.defaultdict(
        collections.Counter
    )
    for record in records:
        by_type[record.task_type][status[record.asset_id]] += 1
    print("== byte-level status by task type (identical / changed / missing)")
    for task_type, counts in sorted(by_type.items()):
        print(
            f"  {task_type:42} {counts['identical']:4} {counts['changed']:4} "
            f"{counts['missing']:4}"
        )
    totals = collections.Counter(status.values())
    print(f"  {'TOTAL':42} {totals['identical']:4} {totals['changed']:4} {totals['missing']:4}")
    missing = sorted(a for a, s in status.items() if s == "missing")
    if missing:
        print(f"missing after rebuild: {missing}")
    print("\n== changed assets (changed files)")
    for asset_id, names in sorted(changed_files.items()):
        shown = ", ".join(names[:6]) + (" ..." if len(names) > 6 else "")
        print(f"  {asset_id}: {shown}")
    return set(changed_files)


def compare_keyed(
    label: str, old: pl.DataFrame, new: pl.DataFrame, key: Sequence[str]
) -> pl.DataFrame:
    """Rows on one side only, and the largest absolute change per shared numeric column."""
    only_old = old.join(new.select(key), on=key, how="anti").height
    only_new = new.join(old.select(key), on=key, how="anti").height
    joined = old.join(new, on=key, how="inner", suffix="_new")
    numeric = [
        c
        for c in old.columns
        if c not in key
        and c + "_new" in joined.columns
        and old.schema[c].is_numeric()
        and new.schema[c].is_numeric()
    ]
    print(
        f"  [{label}] rows old={old.height} new={new.height} shared={joined.height} "
        f"only_old={only_old} only_new={only_new}"
    )
    if numeric and joined.height:
        deltas = joined.select(
            [
                (pl.col(c + "_new").cast(pl.Float64) - pl.col(c).cast(pl.Float64))
                .abs()
                .max()
                .alias(c)
                for c in numeric
            ]
        ).row(0, named=True)
        print(
            "    max |delta|: "
            + ", ".join(f"{c}={v:.3g}" for c, v in deltas.items() if v is not None)
        )
    return joined


def report_significance_flips(
    joined: pl.DataFrame, key: Sequence[str], p_col: str
) -> None:
    if p_col not in joined.columns:
        return
    flipped = joined.filter(
        (pl.col(p_col) < SIGNIFICANCE) != (pl.col(p_col + "_new") < SIGNIFICANCE)
    )
    print(
        f"    significant at {p_col}<{SIGNIFICANCE}: "
        f"old={joined.filter(pl.col(p_col) < SIGNIFICANCE).height} "
        f"new={joined.filter(pl.col(p_col + '_new') < SIGNIFICANCE).height} "
        f"flipped={flipped.height}"
    )
    if flipped.height:
        print(flipped.select(*key, p_col, p_col + "_new"))


def read_magma_out(path: Path) -> pl.DataFrame:
    lines = [line.split() for line in path.read_text().splitlines()]
    rows = [fields for fields in lines if fields and not fields[0].startswith("#")]
    header, body = rows[0], rows[1:]
    frame = pl.DataFrame(body, schema=header, orient="row")
    return frame.with_columns(
        pl.col(c).cast(pl.Float64, strict=False)
        for c in header
        if c not in ("VARIABLE", "TYPE", "FULL_NAME", "GENE")
    )


def compare_csv(record: Record, key: Sequence[str], p_col: str) -> None:
    for name in record.file_md5:
        if name.endswith(".csv"):
            old = pl.read_csv(record.old(name), infer_schema_length=None)
            new = pl.read_csv(record.new(name), infer_schema_length=None)
            joined = compare_keyed(f"{record.asset_id}/{name}", old, new, key)
            report_significance_flips(joined, key, p_col)


def compare_ldsc_h2(record: Record) -> None:
    for name in record.file_md5:
        print(f"  [{record.asset_id}] old:")
        print(pl.read_csv(record.old(name)))
        print("  new:")
        print(pl.read_csv(record.new(name)))


def compare_s_ldsc(record: Record) -> None:
    p_col = (
        CORRECTED_P_COL
        if record.task_type == "MultipleTestingTableTask"
        else "Coefficient_P_value"
    )
    compare_csv(record, key=["Name"], p_col=p_col)


def compare_multiple_testing(record: Record) -> None:
    if "s_ldsc" in record.asset_id:
        compare_s_ldsc(record)
    elif "mr" in record.asset_id:
        compare_csv(record, key=["exposure", "method"], p_col=CORRECTED_P_COL)
    else:
        compare_csv(record, key=["VARIABLE"], p_col=CORRECTED_P_COL)


def compare_mr(record: Record) -> None:
    old = pl.read_csv(record.old("mr_result.csv"), infer_schema_length=None)
    new = pl.read_csv(record.new("mr_result.csv"), infer_schema_length=None)
    key = ["exposure", "method"]
    joined = compare_keyed(f"{record.asset_id}/mr_result.csv", old, new, key)
    report_significance_flips(joined, key, "pval")
    sign = joined.filter(pl.col("b").sign() != pl.col("b_new").sign())
    print(f"    effect sign flips: {sign.height}")
    if sign.height:
        print(sign.select(*key, "b", "b_new", "pval", "pval_new"))


def compare_magma(record: Record) -> None:
    for name in record.file_md5:
        if name.endswith(".gsa.out"):
            key = ["VARIABLE"]
        elif name.endswith(".genes.out"):
            key = ["GENE"]
        else:
            continue
        old = read_magma_out(record.old(name))
        new = read_magma_out(record.new(name))
        # A conditional analysis reports each variable once per model.
        if "MODEL" in old.columns:
            key = ["MODEL", *key]
        joined = compare_keyed(f"{record.asset_id}/{name}", old, new, key)
        report_significance_flips(joined, key, "P")


def credible_sets(path: Path) -> dict[str, set[tuple[object, ...]]]:
    """Variants per credible set; a run that found no credible set writes no file."""
    if not path.exists():
        return {}
    frame = pl.read_parquet(path)
    return {
        str(cs): set(group.select(CS_VARIANT_COLS).iter_rows())
        for (cs,), group in frame.group_by("cs")
    }


def read_sets_yaml(path: Path) -> dict[str, object]:
    return yaml.safe_load(path.read_text()) if path.exists() else {}


def compare_susie(record: Record) -> None:
    old_pip = pl.read_parquet(record.old("pip.parquet"))
    new_pip = pl.read_parquet(record.new("pip.parquet"))
    joined = compare_keyed(record.asset_id + "/pip", old_pip, new_pip, [PIP_ID_COL])
    moved = joined.filter((pl.col("PIP_new") - pl.col("PIP")).abs() > 0.05)
    if moved.height:
        print(f"    variants with |delta PIP| > 0.05: {moved.height}")
        print(moved.sort((pl.col("PIP_new") - pl.col("PIP")).abs(), descending=True).head(5))
    old_meta = read_sets_yaml(record.old("sets.yaml"))
    new_meta = read_sets_yaml(record.new("sets.yaml"))
    old_cs = credible_sets(record.old("combined_cs.parquet"))
    new_cs = credible_sets(record.new("combined_cs.parquet"))
    print(
        f"    credible sets old={sorted(old_cs)} new={sorted(new_cs)}; "
        f"coverage old={old_meta.get('coverage')} new={new_meta.get('coverage')}"
    )
    for cs in sorted(set(old_cs) | set(new_cs)):
        before, after = old_cs.get(cs, set()), new_cs.get(cs, set())
        if before != after:
            print(
                f"    {cs}: size {len(before)} -> {len(after)}, "
                f"shared {len(before & after)}"
            )


def compare_harmonized(record: Record) -> None:
    (name,) = record.file_md5
    columns = [
        GWASLAB_SNPID_COL,
        GWASLAB_EFFECT_ALLELE_COL,
        GWASLAB_BETA_COL,
        GWASLAB_EFFECT_ALLELE_FREQ_COL,
    ]
    old = pl.read_parquet(record.old(name), columns=columns)
    new = pl.read_parquet(record.new(name), columns=columns)
    key = [GWASLAB_SNPID_COL]
    joined = compare_keyed(f"{record.asset_id}", old, new, key)
    swapped = joined.filter(
        pl.col(GWASLAB_EFFECT_ALLELE_COL).cast(pl.String)
        != pl.col(GWASLAB_EFFECT_ALLELE_COL + "_new").cast(pl.String)
    ).height
    sign = joined.filter(
        pl.col(GWASLAB_BETA_COL).sign() != pl.col(GWASLAB_BETA_COL + "_new").sign()
    ).height
    print(f"    effect allele changed: {swapped}; BETA sign changed: {sign}")


def compare_snpvar(record: Record) -> None:
    old = json.loads(record.old("diagnostics.json").read_text())
    new = json.loads(record.new("diagnostics.json").read_text())
    for field in ("alpha_even", "alpha_odd"):
        print(f"  [{record.asset_id}] {field}: {old.get(field)} -> {new.get(field)}")
    old_snpvar = pl.read_parquet(record.old("snpvar.parquet"))
    new_snpvar = pl.read_parquet(record.new("snpvar.parquet"))
    key = [c for c in old_snpvar.columns if not old_snpvar.schema[c].is_float()]
    compare_keyed(record.asset_id + "/snpvar", old_snpvar, new_snpvar, key)


COMPARERS: dict[str, Callable[[Record], None]] = {
    "GenomeReferenceHarmonizationTask": compare_harmonized,
    "SNPHeritabilityByLDSCTask": compare_ldsc_h2,
    "CellAnalysisByLDSCTask": compare_s_ldsc,
    "MultipleTestingTableTask": compare_multiple_testing,
    "TwoSampleMRTask": compare_mr,
    "MagmaGeneSetAnalysisTask": compare_magma,
    "MagmaGeneAnalysisTask": compare_magma,
    "SusieRFinemapTask": compare_susie,
    "L2RegularizedSldscSnpvarTask": compare_snpvar,
}


def main() -> None:
    records = load_records()
    changed = report_md5_status(records)
    with pl.Config(tbl_rows=20, tbl_cols=12, tbl_width_chars=200):
        for task_type, comparer in COMPARERS.items():
            selected = [
                r
                for r in records
                if r.task_type == task_type and r.asset_id in changed and r.copied
            ]
            if not selected:
                continue
            print(f"\n== {task_type}: content of {len(selected)} changed assets")
            for record in sorted(selected, key=lambda r: r.asset_id):
                comparer(record)


if __name__ == "__main__":
    main()
