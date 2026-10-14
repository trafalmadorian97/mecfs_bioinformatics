"""
Confirm that the zero-SE variants in the Kerrebijn fibromyalgia sumstats are what
break cross-trait LDSC, by running the Fibromyalgia x Schizophrenia rg twice:

  run A: data exactly as the CT-LDSC task feeds it today  -> expected to fail
  run B: same data with SE == 0 variants dropped          -> expected to succeed

gwaslab catches the LinAlgError internally and then crashes in its own error
handler (traceback.format_exc(ex) passes the exception as `limit`), so run A
surfaces as a TypeError rather than the LinAlgError itself.

Run:
    pixi r python experiments/claude/fibro_ct_ldsc_zero_se_repro.py \
        2>&1 | tee experiments/claude/logs/fibro_ct_ldsc_zero_se_repro.log
"""

import copy
from pathlib import Path

import gwaslab

from mecfs_bio.build_system.task.gwaslab.gwaslab_genetic_corr_by_ct_ldsc_task import (
    FilterSettings,
    filter_sumstats,
    get_compatible_snps_polars,
)
from mecfs_bio.constants.gwaslab_constants import GWASLAB_RSID_COL

ASSET_STORE = Path("assets/base_asset_store")
FIBRO_PICKLE = (
    ASSET_STORE
    / "gwas/fibromyalgia/kerrebijn_et_al/gwaslab_sumstats/kerrebijin_fibro_sumstats_37.pickle"
)
SCZ_PICKLE = (
    ASSET_STORE
    / "gwas/schizophrenia/pgc_2022/gwaslab_sumstats/pgc_2022_sch_sumstats_37.pickle"
)
LD_DIR = (
    ASSET_STORE
    / "reference_data/linkage_disequilibrium_scores/thousand_genomes_phase_3_v1/extracted"
    / "thousand_genomes_phase_3_v1_eur_ld_scores_extracted"
)

RG_KWARGS = dict(
    ref_ld_chr=str(LD_DIR) + "/LDscore.@",
    w_ld_chr=str(LD_DIR) + "/LDscore.@",
    build="19",
    samp_prev="0.5,0.408",
    pop_prev="0.027,0.01",
)


def load_filtered(path: Path) -> gwaslab.Sumstats:
    sumstats = gwaslab.load_pickle(str(path))
    sumstats.infer_build()
    filter_sumstats(sumstats, FilterSettings(), build="19")
    return sumstats


def run_rg(fibro: gwaslab.Sumstats, scz: gwaslab.Sumstats, label: str) -> None:
    print(f"\n############ {label} ############")
    print(f"fibro rows={len(fibro.data)}  zero_SE={int((fibro.data['SE'] == 0).sum())}")

    compatible = get_compatible_snps_polars(fibro.data, scz.data)
    scz = copy.deepcopy(scz)
    scz.data = scz.data.loc[
        scz.data[GWASLAB_RSID_COL].isin(compatible[GWASLAB_RSID_COL])
    ]

    try:
        fibro.estimate_rg_by_ldsc(other_traits=[scz], rg="Fibromyalgia,Schizophrenia", **RG_KWARGS)
        print(f"RESULT [{label}]: SUCCESS")
        print(fibro.ldsc_rg.to_string())
    except Exception as exc:  # noqa: BLE001 - we are characterising the failure
        print(f"RESULT [{label}]: FAILED with {type(exc).__name__}: {exc}")


def main() -> None:
    fibro = load_filtered(FIBRO_PICKLE)
    scz = load_filtered(SCZ_PICKLE)

    run_rg(copy.deepcopy(fibro), scz, "A: as-is (SE == 0 kept)")

    fibro_clean = copy.deepcopy(fibro)
    fibro_clean.data = fibro_clean.data.loc[fibro_clean.data["SE"] > 0]
    run_rg(fibro_clean, scz, "B: SE == 0 dropped")


if __name__ == "__main__":
    main()
