"""
Compare DecodeME build-37 downstream results before vs after switching to genome-reference
harmonization. "before" is a snapshot taken from the asset store prior to the rebuild;
"after" is read live from the asset store once rerun_downstream.py has finished.

Reports, per analysis, only what matters qualitatively:
- S-LDSC: number of significant (rejected-null) cell types/tissues per panel, and the
  top-5 by coefficient p-value, before vs after.
- MAGMA GTEx / HBA: number of Bonferroni-significant sets and the top-5 by p-value.
- LDSC: liability-scale h2, its SE, intercept, ratio.
- SUSIE (non-PolyFun, with-palindromes loci): the top-PIP variant of each credible set.

Run (after rerun_downstream.py completes):
  pixi r python -m experiments.claude.genome_reference_harmonization.compare_downstream \
    --before <snapshot_dir> \
    2>&1 | tee experiments/claude/genome_reference_harmonization/compare_downstream.log
"""

import argparse
from pathlib import Path

import polars as pl

STORE = Path("assets/base_asset_store/gwas/ME_CFS/DecodeME")
AFTER_ANALYSIS = STORE / "analysis"
AFTER_MAGMA = STORE / "analysis_results" / "magma"

# The six aggregate build-37 S-LDSC panels (excludes build_38 and the c0..c8 experiments).
SLDSC_PANELS = [
    "decode_me_gwas_1_cahoy_cns_s_ldsc_multiple_testing_correction.csv",
    "decode_me_gwas_1_corces_atac_s_ldsc_multiple_testing_correction.csv",
    "decode_me_gwas_1_gtex_brain_s_ldsc_multiple_testing_correction.csv",
    "decode_me_gwas_1_immgen_s_ldsc_multiple_testing_correction.csv",
    "decode_me_gwas_1_multi_tissue_chromatin_s_ldsc_multiple_testing_correction.csv",
    "decode_me_gwas_1_multi_tissue_gene_expression_s_ldsc_multiple_testing_correction.csv",
]

MAGMA_GSA = {
    "MAGMA GTEx specific tissue": (
        "decode_me_gwas_1_build_37_magma_ensemble_specific_tissue_gene_covar_analysis"
        "/gene_set_analysis_output.gsa.out"
    ),
    "MAGMA HBA gene covar": (
        "decode_me_hba_magma_tasks_hba_gene_covar/gene_set_analysis_output.gsa.out"
    ),
    "MAGMA HBA conditional": (
        "decode_me_hba_magma_tasks_hba_magma_conditional_analysis"
        "/gene_set_analysis_output.gsa.out"
    ),
}


def _fmt(x: float) -> str:
    return f"{x:.3g}"


def compare_sldsc(before: Path) -> None:
    print("\n" + "=" * 78)
    print("S-LDSC: significant cell types/tissues (Reject Null) and top-5 by p-value")
    print("=" * 78)
    for panel in SLDSC_PANELS:
        b_path, a_path = before / "analysis" / panel, AFTER_ANALYSIS / panel
        name = panel.replace("decode_me_gwas_1_", "").replace(
            "_s_ldsc_multiple_testing_correction.csv", ""
        )
        if not (b_path.exists() and a_path.exists()):
            print(f"\n[{name}] MISSING before={b_path.exists()} after={a_path.exists()}")
            continue
        b, a = pl.read_csv(b_path), pl.read_csv(a_path)
        for label, df in (("before", b), ("after", a)):
            sig = df.filter(pl.col("Reject Null")).height
            top = (
                df.sort("Coefficient_P_value")
                .head(5)
                .select("Name", "Coefficient_P_value", "_Corrected P Value_")
            )
            print(f"\n[{name}] {label}: {sig} significant / {df.height} tested")
            for row in top.iter_rows():
                print(f"    {row[0][:48]:48s} p={_fmt(row[1])}  q={_fmt(row[2])}")


def _read_gsa(path: Path) -> pl.DataFrame:
    lines = [
        ln
        for ln in path.read_text().splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    header = lines[0].split()
    rows = [ln.split(maxsplit=len(header) - 1) for ln in lines[1:]]
    df = pl.DataFrame(rows, schema=header, orient="row")
    return df.with_columns(pl.col("P").cast(pl.Float64))


def compare_magma(before: Path) -> None:
    print("\n" + "=" * 78)
    print("MAGMA: Bonferroni-significant sets and top-5 by p-value")
    print("=" * 78)
    for name, rel in MAGMA_GSA.items():
        b_path, a_path = before / "magma" / rel, AFTER_MAGMA / rel
        if not (b_path.exists() and a_path.exists()):
            print(f"\n[{name}] MISSING before={b_path.exists()} after={a_path.exists()}")
            continue
        b, a = _read_gsa(b_path), _read_gsa(a_path)
        for label, df in (("before", b), ("after", a)):
            n = df.height
            thresh = 0.05 / n
            sig = df.filter(pl.col("P") < thresh).height
            top = df.sort("P").head(5).select("VARIABLE", "P")
            print(f"\n[{name}] {label}: {sig} Bonferroni-sig / {n} sets (a=0.05/{n})")
            for var, p in top.iter_rows():
                print(f"    {var[:52]:52s} p={_fmt(p)}")


def compare_ldsc(before: Path) -> None:
    print("\n" + "=" * 78)
    print("LDSC: liability-scale SNP heritability")
    print("=" * 78)
    fn = "decode_me_gwas_1_heritability_by_ldsc.csv"
    b_path, a_path = before / "analysis" / fn, AFTER_ANALYSIS / fn
    if not (b_path.exists() and a_path.exists()):
        print(f"MISSING before={b_path.exists()} after={a_path.exists()}")
        return
    cols = ["h2_liab", "h2_liab_se", "Lambda_gc", "Mean_chi2", "Intercept", "Ratio"]
    b, a = pl.read_csv(b_path), pl.read_csv(a_path)
    for c in cols:
        bv, av = b[c][0], a[c][0]
        print(f"    {c:14s} before={bv!s:14s} after={av!s:14s}")


def _top_variant_per_cs(path: Path) -> dict[str, tuple[str, float]]:
    df = pl.read_parquet(path)
    out: dict[str, tuple[str, float]] = {}
    for cs, sub in df.group_by("cs"):
        top = sub.sort("PIP", descending=True).row(0, named=True)
        key = cs[0] if isinstance(cs, tuple) else cs
        out[str(key)] = (
            f"{top['CHR']}:{top['POS']}:{top['EA']}:{top['NEA']}",
            float(top["PIP"]),
        )
    return out


def compare_susie(before: Path) -> None:
    print("\n" + "=" * 78)
    print("SUSIE (non-PolyFun, with-palindromes loci): top-PIP variant per credible set")
    print("=" * 78)
    pattern = "*palindromes_keep*cs_from_directory.parquet"
    names = sorted(
        {p.name for p in (before / "analysis").glob(pattern)}
        | {p.name for p in AFTER_ANALYSIS.glob(pattern)}
    )
    for fn in names:
        b_path, a_path = before / "analysis" / fn, AFTER_ANALYSIS / fn
        short = fn.replace("_copy_cs_from_directory.parquet", "")
        if not (b_path.exists() and a_path.exists()):
            print(f"\n[{short}] MISSING before={b_path.exists()} after={a_path.exists()}")
            continue
        b, a = _top_variant_per_cs(b_path), _top_variant_per_cs(a_path)
        css = sorted(set(b) | set(a))
        changed = any(b.get(c) != a.get(c) for c in css)
        flag = "  <-- CHANGED" if changed else ""
        print(f"\n[{short}]{flag}")
        for c in css:
            bv, av = b.get(c), a.get(c)
            mark = "" if bv == av else "   *"
            print(f"    {c}: before={bv}  after={av}{mark}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--before", type=Path, required=True)
    args = parser.parse_args()
    compare_ldsc(args.before)
    compare_sldsc(args.before)
    compare_magma(args.before)
    compare_susie(args.before)


if __name__ == "__main__":
    main()
