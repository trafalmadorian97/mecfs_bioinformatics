"""
How sensitive is the (now linear) GWAS-by-subtraction to the DecodeME sample
size we pass? Theory says the per-SNP remainder (est, se_c, Z, N_eff) is exactly
invariant to N_D (the N cancels through sumstats -> ldsc -> kernel), with only a
weak 2nd-order effect via the ldsc chi^2-outlier filter (max(0.001*N, 80)).

Rebuild the subtraction at N_D in {58789 (effective N), 275488 (total, baseline),
2754880 (10x total)} and compare the remainder per-SNP outputs and the internal
ldsc S matrix. Cleans up the temp asset dirs afterward.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path

import numpy as np
import polars as pl

from mecfs_bio.analysis.runner.default_runner import DEFAULT_RUNNER
from mecfs_bio.assets.gwas.me_cfs.decode_me.processed_gwas_data.decode_me_annovar_37_rsids_assignment import (
    DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS,
)
from mecfs_bio.assets.gwas.multisite_pain.johnston_et_al.analysis.johnston_standard_analysis import (
    JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS,
)
from mecfs_bio.assets.reference_data.genomic_sem_reference.genomes_1k import (
    GENOMES1K_REFERENCE_FOR_GENOMIC_SEM,
)
from mecfs_bio.assets.reference_data.genomic_sem_reference.hapmap3_snpist import (
    HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
)
from mecfs_bio.assets.reference_data.linkage_disequilibrium_score_reference_data.extracted.eur_ld_scores_thousand_genome_phase_3_v1_extracted import (
    THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
)
from mecfs_bio.build_system.task.pipes.composite_pipe import CompositePipe
from mecfs_bio.build_system.task.pipes.compute_p_pipe import ComputePPipe
from mecfs_bio.build_system.task.pipes.rename_col_pipe import RenameColPipe
from mecfs_bio.build_system.task.r_tasks.genomic_sem._genomic_sem_config import (
    GenomicSEMSumstatsSource,
)
from mecfs_bio.build_system.task.r_tasks.genomic_sem.genomic_sem_gwas_by_subtraction_full_python_task import (
    GenomicSEMGWASBySubtractionFullPythonTask,
    GWASBySubtractionFullPythonConfig,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

ANALYSIS_DIR = Path(
    "assets/base_asset_store/gwas/multi_trait/genomic_sem/analysis"
)
BASELINE = ANALYSIS_DIR / "decode_me_minus_pain_subtraction_full_python_ols"
PAIN_N = 387649

# N_D label -> asset_id suffix.
VARIANTS = {
    "neff_58789": 58789,        # effective N: 4/(1/15579 + 1/259909)
    "tenx_2754880": 2754880,    # 10x the total N
}


def make_task(n_d: int, asset_id: str) -> GenomicSEMGWASBySubtractionFullPythonTask:
    return GenomicSEMGWASBySubtractionFullPythonTask.create(
        asset_id=asset_id,
        composite_trait_source=GenomicSEMSumstatsSource(
            task=DECODE_ME_GWAS_1_37_ANNOVAR_DBSNP150_RSID_ASSIGNED_KEEP_AMBIGUOUS.join_task,
            alias="DECODE_ME",
            sample_size=n_d,
            pipe=CompositePipe(
                [
                    RenameColPipe(old_name="rsid", new_name=GWASLAB_RSID_COL),
                    ComputePPipe(),
                ]
            ),
        ),
        reference_trait_source=GenomicSEMSumstatsSource(
            task=JOHNSTON_ET_AL_PAIN_STANDARD_ANALYSIS.magma_tasks.parquet_file_task,
            alias="Multsite_Pain",
            sample_size=PAIN_N,
        ),
        ld_ref_task=THOUSAND_GENOME_EUR_LD_REFERENCE_DATA_V1_EXTRACTED,
        hapmap_snps_task=HAPMAP3_SNPLIST_FOR_GENOMIC_SEM,
        sumstats_ref_task=GENOMES1K_REFERENCE_FOR_GENOMIC_SEM,
        config=GWASBySubtractionFullPythonConfig(ld_file_filename_pattern="LDscore."),
    )


def _remainder(d: Path) -> pl.DataFrame:
    return pl.read_parquet(d / "gwas_results" / "remainder_factor.parquet").select(
        "SNP", "est", "se_c", "Z_Estimate", "N_eff"
    )


def main() -> None:
    cols = ["est", "se_c", "Z_Estimate", "N_eff"]
    asset_ids = {
        k: f"decode_me_minus_pain_subtraction_full_python_sens_{k}" for k in VARIANTS
    }
    try:
        for k, n_d in VARIANTS.items():
            print(f"=== building N_D={n_d} ({k}) ===", flush=True)
            DEFAULT_RUNNER.run([make_task(n_d, asset_ids[k])])

        base = _remainder(BASELINE)
        print("\n===== internal ldsc S matrix vs N_D (S11 ~ 1/N_D, S12 ~ 1/sqrt(N_D)) =====")
        base_S = np.loadtxt(BASELINE / "ldsc_S.csv", delimiter=",", skiprows=1)
        print(f"N_D=275488 (baseline): S11={base_S[0,0]:.5f} S12={base_S[0,1]:.5f} S22={base_S[1,1]:.5f}")
        for k, n_d in VARIANTS.items():
            S = np.loadtxt(ANALYSIS_DIR / asset_ids[k] / "ldsc_S.csv", delimiter=",", skiprows=1)
            print(f"N_D={n_d:>8} ({k}): S11={S[0,0]:.5f} S12={S[0,1]:.5f} S22={S[1,1]:.5f}  "
                  f"[S11*N_D/M0={S[0,0]*n_d/275488:.5f} vs baseline S11]")

        print("\n===== per-SNP remainder: baseline (N_D=275488) vs variants =====")
        for k, n_d in VARIANTS.items():
            v = _remainder(ANALYSIS_DIR / asset_ids[k])
            j = base.join(v, on="SNP", suffix="_v")
            print(f"\n--- N_D={n_d} ({k}); {j.height} shared SNPs ---")
            for c in cols:
                a = j[c].to_numpy()
                b = j[f"{c}_v"].to_numpy()
                denom = np.maximum(np.abs(a), 1e-12)
                rel = np.abs(a - b) / denom
                print(f"  {c:<11} max|rel|={rel.max():.2e}  median|rel|={np.median(rel):.2e}  "
                      f"corr={np.corrcoef(a, b)[0, 1]:.6f}")
    finally:
        for k in VARIANTS:
            d = ANALYSIS_DIR / asset_ids[k]
            if d.exists():
                shutil.rmtree(d)
                print(f"cleaned up {d}")


if __name__ == "__main__":
    main()
    print("\nDONE", file=sys.stderr)
