"""Summarize the Genebass and FinnGen PheWAS lookups written by phewas_lookups.py.

Multiple-testing thresholds are Bonferroni over the phenotypes tested in each resource for each test
(one gene-level test per burden set in Genebass, one variant in FinnGen).

Run: pixi r python experiments/claude/net_a457p_gnomad/summarize_lookups.py
"""

from pathlib import Path

import polars as pl

OUT_DIR = Path(__file__).parent
ALPHA = 0.05
N_TOP = 8


def summarize_genebass(genebass: pl.DataFrame) -> None:
    print("== Genebass SLC6A2 gene-level PheWAS (UKB 394k exomes; SKAT-O p)")
    for (annotation,), group in genebass.group_by("annotation", maintain_order=True):
        tested = group.filter(pl.col("Pvalue").is_not_null())
        threshold = ALPHA / tested.height
        n_variants = tested["total_variants"].max()
        print(
            f"\n-- {annotation}: phenotypes={tested.height} max_variants={n_variants} Bonferroni={threshold:.1e} "
            f"n_sig={tested.filter(pl.col('Pvalue') < threshold).height} n_p<1e-3={tested.filter(pl.col('Pvalue') < 1e-3).height} "
            f"(expected {tested.height * 1e-3:.1f})"
        )
        for row in tested.sort("Pvalue").head(N_TOP).iter_rows(named=True):
            print(
                f"   p={row['Pvalue']:.2e} burden_p={row['Pvalue_Burden']} beta={row['BETA_Burden']} "
                f"{row['trait_type']}:{row['phenocode']} {str(row['description'])[:60]} {str(row['coding_description'] or '')[:30]} "
                f"cases={row['n_cases']}"
            )


def summarize_finngen(finngen: pl.DataFrame) -> None:
    print(
        "\n== FinnGen R12 single-variant PheWAS (the variant API returns statistics only for nominally associated "
        "endpoints and null otherwise, so Bonferroni and expected counts use all endpoints)"
    )
    for (variant_id, rsid, role), all_endpoints in finngen.group_by(
        "variant_id", "rsid", "role", maintain_order=True
    ):
        group = all_endpoints.filter(pl.col("pval").is_not_null())
        threshold = ALPHA / all_endpoints.height
        print(
            f"\n-- {rsid} {variant_id} [{role}] FinnGen alt AF={group['maf'].median():.4f} "
            f"endpoints={all_endpoints.height} reported={group.height} "
            f"Bonferroni={threshold:.1e} n_sig={group.filter(pl.col('pval') < threshold).height} "
            f"n_p<5e-8={group.filter(pl.col('pval') < 5e-8).height} n_p<1e-3={group.filter(pl.col('pval') < 1e-3).height} "
            f"(expected {all_endpoints.height * 1e-3:.1f})"
        )
        for row in group.sort("pval").head(N_TOP).iter_rows(named=True):
            print(
                f"   p={row['pval']:.2e} beta={row['beta']:.3f} se={row['sebeta']:.3f} {row['phenocode']} "
                f"{row['phenostring'][:60]} cases={row['n_case']}"
            )


def summarize_ot(ot_sets: pl.DataFrame) -> None:
    print("\n== Open Targets GWAS credible sets containing a lookup variant")
    for row in ot_sets.sort("lead_p").iter_rows(named=True):
        print(
            f"   {row['rsid']} pip={row['lookup_variant_pip']} trait={row['trait'][:60]} {row['first_author']} "
            f"N={row['n_samples']} lead={row['lead_variant']} p={row['lead_p']:.1e} L2G_top={row['l2g_top_gene']}"
        )


def main() -> None:
    summarize_genebass(pl.read_parquet(OUT_DIR / "genebass_slc6a2.parquet"))
    summarize_finngen(pl.read_parquet(OUT_DIR / "finngen_slc6a2_lookups.parquet"))
    summarize_ot(
        pl.read_parquet(OUT_DIR / "ot_gwas_credible_sets_containing_lookups.parquet")
    )


if __name__ == "__main__":
    main()
