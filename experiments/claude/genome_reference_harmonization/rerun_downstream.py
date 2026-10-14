"""
Rerun the main build-37 DecodeME downstream analyses that depend on the harmonization
output, so we can check whether switching to genome-reference harmonization changed any
qualitative result.

Covered (all consume the switched annovar-37 rsID-assignment chains):
- MAGMA GTEx specific-tissue bar plot + its gene-set analysis
- MAGMA human-brain-atlas (HBA) analysis
- S-LDSC (all terminal tables)
- univariate LDSC SNP heritability
- non-PolyFun SUSIE fine-mapping (the six with-palindromes loci)

PolyFun SUSIE is intentionally skipped (known outstanding fixes).

All three entrypoints share one on-disk verifying-trace cache, so the common upstream
(harmonization -> rsID assignment) is built once and reused. Rebuilds are driven by the
trace: the production *_genome_reference_harmonized assets do not exist yet, so the
cascade fires without must_rebuild_transitive.

Run:
  pixi r python -m experiments.claude.genome_reference_harmonization.rerun_downstream \
    2>&1 | tee experiments/claude/genome_reference_harmonization/rerun_downstream.log
"""

from mecfs_bio.analysis.decode_me_fine_mapping import (
    run_fine_mapping_decode_me_analysis,
)
from mecfs_bio.analysis.decode_me_initial_analysis import (
    run_initial_decode_me_analysis,
)
from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.analysis.decode_me_gwas_1_ldsc import (
    DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD,
)


def main() -> None:
    print("=== [1/3] initial analysis: MAGMA GTEx + MAGMA HBA + S-LDSC ===", flush=True)
    run_initial_decode_me_analysis()

    print("=== [2/3] univariate LDSC SNP heritability ===", flush=True)
    DEFAULT_RUNNER.run(
        [DECODE_ME_GWAS_1_HERITABILITY_BY_LDSC_MD],
        incremental_save=True,
        must_rebuild_transitive=[],
    )

    print("=== [3/3] non-PolyFun SUSIE fine-mapping (with-palindromes loci) ===", flush=True)
    run_fine_mapping_decode_me_analysis()

    print("=== done ===", flush=True)


if __name__ == "__main__":
    main()
